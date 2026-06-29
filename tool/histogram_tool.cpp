#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include "ExportReader.h"
#include "FinalDrainReader.h"
#include "FlatConcurrentHashMap.h"
#include "HighFrequencyInsertThreadPool.h"
#include "LowFrequencyQueryThreadPool.h"

#include "../src/RingMemoryPool.h"

namespace
{
    constexpr uint64_t HISTOGRAM_RING_CAPACITY = 1ULL << 7;
    constexpr uint64_t HISTOGRAM_BLOCK_BYTES = 128ULL * 1024ULL;
    constexpr uint64_t ROOT_BUCKET_COUNT = 1ULL << (2 * ROOT_BASES);
    constexpr uint64_t BYTES_PER_GIB = 1024ULL * 1024ULL * 1024ULL;

    struct Options
    {
        std::string tmp_dir;
        uint32_t k_len = 0;
        uint32_t max_threads = 0;
        uint64_t max_memory_bytes = 0;
        std::string output_file;
        uint32_t min_freq = 1;
        uint32_t max_freq = 10000;
    };

    struct RootFileInfo
    {
        std::string filename;
        uint64_t record_count = 0;
        uint64_t file_size = 0;
    };

    using AtomicHistogram = std::vector<std::atomic<int64_t>>;

    class ProgressPrinter
    {
    public:
        explicit ProgressPrinter(const uint64_t total_bytes)
            : total_bytes_(total_bytes)
        {
        }

        void start()
        {
            std::cout << "progress 0%";
            std::flush(std::cout);
            if (total_bytes_ == 0)
            {
                last_percent_ = 100;
            }
        }

        void add(const uint64_t bytes)
        {
            if (total_bytes_ == 0 || bytes == 0)
            {
                return;
            }

            completed_bytes_ += bytes;
            uint64_t percent = (completed_bytes_ >= total_bytes_)
                ? 100
                : (completed_bytes_ * 100) / total_bytes_;
            if (percent > 99)
            {
                percent = 99;
            }

            if (percent > last_percent_)
            {
                std::cout << "\rprogress " << percent << "%";
                std::flush(std::cout);
            }
            last_percent_ = percent;
        }

        void finish()
        {
            if (!finished_)
            {
                std::cout << "\rprogress 100%" << std::endl;
                finished_ = true;
            }
        }

    private:
        uint64_t total_bytes_ = 0;
        uint64_t completed_bytes_ = 0;
        uint64_t last_percent_ = 0;
        bool finished_ = false;
    };

    void init_histogram(AtomicHistogram& histogram)
    {
        for (auto& value : histogram)
        {
            value.store(0, std::memory_order_relaxed);
        }
    }

    [[nodiscard]] std::string with_trailing_slash(std::string path)
    {
        if (path.empty())
        {
            return "./";
        }
        const char tail = path.back();
        if (tail != '/' && tail != '\\')
        {
            path.push_back('/');
        }
        return path;
    }

    [[nodiscard]] uint32_t parse_u32(const char* text, const char* name)
    {
        size_t parsed = 0;
        const unsigned long long value = std::stoull(text, &parsed);
        if (parsed != std::strlen(text) || value > std::numeric_limits<uint32_t>::max()) [[unlikely]]
        {
            std::cerr << "invalid " << name << ": " << text << std::endl;
            exit(-1);
        }
        return static_cast<uint32_t>(value);
    }

    [[nodiscard]] uint64_t parse_memory_gib(const char* text)
    {
        size_t parsed = 0;
        const long double value = std::stold(text, &parsed);
        if (parsed != std::strlen(text) || value <= 0.0L) [[unlikely]]
        {
            std::cerr << "invalid max_memory_gb: " << text << std::endl;
            exit(-1);
        }

        const long double bytes = value * static_cast<long double>(BYTES_PER_GIB);
        if (bytes > static_cast<long double>(std::numeric_limits<uint64_t>::max())) [[unlikely]]
        {
            std::cerr << "max_memory_gb is too large" << std::endl;
            exit(-1);
        }
        return static_cast<uint64_t>(bytes);
    }

    [[nodiscard]] Options parse_options(int argc, char* argv[])
    {
        if (argc < 6 || argc > 8)
        {
            std::cerr << "Usage: histogram_tool <tmp_dir> <k_len> <max_threads> <max_memory_gb> <output_file> [min_freq=1] [max_freq=10000]"
                << std::endl;
            exit(-1);
        }

        Options options;
        options.tmp_dir = with_trailing_slash(argv[1]);
        options.k_len = parse_u32(argv[2], "k_len");
        options.max_threads = parse_u32(argv[3], "max_threads");
        options.max_memory_bytes = parse_memory_gib(argv[4]);
        options.output_file = argv[5];
        if (argc >= 7)
        {
            options.min_freq = parse_u32(argv[6], "min_freq");
        }
        if (argc >= 8)
        {
            options.max_freq = parse_u32(argv[7], "max_freq");
        }

        if (options.k_len == 0 || options.k_len > MAX_K)
        {
            std::cerr << "invalid k_len: must be in [1, MAX_K]" << std::endl;
            exit(-1);
        }
        if (options.max_threads == 0)
        {
            std::cerr << "invalid max_threads: must be greater than 0" << std::endl;
            exit(-1);
        }
        if (options.max_freq < options.min_freq)
        {
            std::cerr << "invalid max_freq: must be greater than or equal to min_freq" << std::endl;
            exit(-1);
        }
        return options;
    }

    [[nodiscard]] std::string root_filename(const std::string& tmp_dir, const uint64_t root_id)
    {
        return tmp_dir + "root_" + std::to_string(root_id) + ".bin";
    }

    void close_fd(const int fd) noexcept
    {
        if (fd >= 0)
        {
            ::close(fd);
        }
    }

    template <uint32_t N>
    [[nodiscard]] std::vector<RootFileInfo> collect_root_files(const std::string& tmp_dir, uint64_t& expected_unique_insert)
    {
        std::vector<RootFileInfo> root_files;
        root_files.reserve(ROOT_BUCKET_COUNT);
        expected_unique_insert = 0;

        for (uint64_t root_id = 0; root_id < ROOT_BUCKET_COUNT; ++root_id)
        {
            const std::string filename = root_filename(tmp_dir, root_id);
            const int fd = ::open(filename.c_str(), O_RDONLY);
            if (fd < 0)
            {
                if (errno == ENOENT)
                {
                    continue;
                }
                std::cerr << "failed to open " << filename << ": " << std::strerror(errno) << std::endl;
                exit(-1);
            }

            struct stat st
            {
            };
            if (::fstat(fd, &st) != 0)
            {
                close_fd(fd);
                std::cerr << "failed to stat " << filename << ": " << std::strerror(errno) << std::endl;
                exit(-1);
            }
            close_fd(fd);

            const uint64_t file_size = static_cast<uint64_t>(st.st_size);
            if (file_size % sizeof(ExportRecord<N>) != 0)
            {
                std::cerr << "invalid root file size for " << filename << std::endl;
                exit(-1);
            }

            const uint64_t record_count = file_size / sizeof(ExportRecord<N>);
            expected_unique_insert = expected_unique_insert + record_count;
            if (record_count > 0)
            {
                root_files.push_back(RootFileInfo{ filename, record_count, file_size });
            }
        }

        return root_files;
    }

    template <uint32_t N>
    [[nodiscard]] uint64_t packed_kmer_bytes_for_k(const uint32_t k_len)
    {
        const uint64_t full_data_count = k_len / BASES_PER_U64T;
        const uint64_t tail_bits = 2ULL * (k_len % BASES_PER_U64T);
        const uint64_t tail_bytes = (tail_bits + 7ULL) / 8ULL;
        const uint64_t packed_kmer_bytes = full_data_count * sizeof(uint64_t) + tail_bytes;
        if (packed_kmer_bytes == 0 || k_len > N * BASES_PER_U64T) [[unlikely]]
        {
            std::cerr << "invalid k-mer length for packed byte calculation: " << k_len << std::endl;
            exit(-1);
        }
        return packed_kmer_bytes;
    }

    template <uint32_t N>
    [[nodiscard]] uint64_t low_file_size_bytes(const std::string& tmp_dir, const uint32_t k_len)
    {
        const std::string filename = tmp_dir + "low.bin";
        const int fd = ::open(filename.c_str(), O_RDONLY);
        if (fd < 0) [[unlikely]]
        {
            std::cerr << "failed to open " << filename << ": " << std::strerror(errno) << std::endl;
            exit(-1);
        }

        struct stat st
        {
        };
        if (::fstat(fd, &st) != 0) [[unlikely]]
        {
            close_fd(fd);
            std::cerr << "failed to stat " << filename << ": " << std::strerror(errno) << std::endl;
            exit(-1);
        }
        close_fd(fd);

        const uint64_t file_size = static_cast<uint64_t>(st.st_size);
        const uint64_t packed_kmer_bytes = packed_kmer_bytes_for_k<N>(k_len);
        if (file_size % packed_kmer_bytes != 0) [[unlikely]]
        {
            std::cerr << "invalid low.bin size: " << file_size
                << " is not divisible by packed k-mer size "
                << packed_kmer_bytes << std::endl;
            exit(-1);
        }
        return file_size;
    }

    template <uint32_t N>
    [[nodiscard]] uint64_t estimate_peak_memory_bytes(
        const uint64_t expected_unique_insert,
        const uint32_t max_threads,
        const size_t hist_size)
    {
        uint64_t total = FlatConcurrentHashMap<N>::required_mmap_bytes(expected_unique_insert);

        total += HISTOGRAM_RING_CAPACITY * HISTOGRAM_BLOCK_BYTES;
        total += HISTOGRAM_RING_CAPACITY * HISTOGRAM_BLOCK_BYTES;

        total += static_cast<uint64_t>(hist_size) * sizeof(std::atomic<int64_t>);
        total += static_cast<uint64_t>(hist_size) * sizeof(int64_t) * max_threads;

        constexpr uint64_t io_slack_bytes = 8ULL * 1024ULL * 1024ULL;
        total += io_slack_bytes;

        return total;
    }

    template <uint32_t N>
    void enqueue_high_records(
        const std::vector<RootFileInfo>& root_files,
        SPMCRingMemoryPool<HISTOGRAM_RING_CAPACITY>& pool,
        ProgressPrinter* progress)
    {
        using Record = ExportRecord<N>;
        constexpr uint64_t RECORDS_PER_BLOCK = HISTOGRAM_BLOCK_BYTES / sizeof(Record);
        static_assert(RECORDS_PER_BLOCK > 0, "histogram high-frequency block is too small");

        for (const RootFileInfo& root_file : root_files)
        {
            FinalDrainReader<N> reader;
            reader.open(root_file.filename);

            while (!reader.finished())
            {
                char* block_ptr = nullptr;
                pool.producer_dequeue(block_ptr);
                auto* records = reinterpret_cast<Record*>(block_ptr);
                const uint64_t count = reader.read_records(records, RECORDS_PER_BLOCK);
                if (count == 0) [[unlikely]]
                {
                    pool.consumer_enqueue(block_ptr);
                    break;
                }
                progress->add(count * sizeof(Record));

                pool.producer_enqueue(content_type{ block_ptr, count });
            }
        }

        pool.producer_set_finished();
    }

    template <uint32_t N>
    void enqueue_low_kmers(
        const uint32_t k_len,
        SPMCRingMemoryPool<HISTOGRAM_RING_CAPACITY>& pool,
        const uint64_t packed_kmer_bytes,
        ProgressPrinter* progress)
    {
        const uint64_t kmers_per_block = HISTOGRAM_BLOCK_BYTES / packed_kmer_bytes;
        if (kmers_per_block == 0) [[unlikely]]
        {
            std::cerr << "histogram low-frequency block is too small for packed k-mer size "
                << packed_kmer_bytes << std::endl;
            exit(-1);
        }

        ExportReader<N> reader(k_len);
        reader.open(temp_dir + "low.bin");

        while (!reader.finished())
        {
            char* block_ptr = nullptr;
            pool.producer_dequeue(block_ptr);
            const uint64_t count = reader.read_packed_kmers(block_ptr, kmers_per_block);
            if (count == 0) [[unlikely]]
            {
                pool.consumer_enqueue(block_ptr);
                break;
            }
            progress->add(count * packed_kmer_bytes);

            pool.producer_enqueue(content_type{ block_ptr, count });
        }

        pool.producer_set_finished();
    }

    void write_histogram(
        const std::string& output_filename,
        const AtomicHistogram& histogram,
        const uint32_t min_freq)
    {
        std::ofstream output(output_filename, std::ios::out | std::ios::trunc);
        if (!output)
        {
            std::cerr << "failed to open output file: " << output_filename << std::endl;
            exit(1);
        }

        for (size_t i = 0; i < histogram.size(); ++i)
        {
            output << (static_cast<uint64_t>(min_freq) + i) << '\t'
                << histogram[i].load(std::memory_order_relaxed) << '\n';
        }
    }

    template <uint32_t N>
    int run_histogram_tool(const Options& options)
    {
        temp_dir = options.tmp_dir;

        const size_t hist_size = static_cast<size_t>(
            static_cast<uint64_t>(options.max_freq) - static_cast<uint64_t>(options.min_freq) + 1ULL);
        const uint32_t worker_count = std::max<uint32_t>(1, options.max_threads);

        uint64_t expected_unique_insert = 0;
        const std::vector<RootFileInfo> root_files = collect_root_files<N>(options.tmp_dir, expected_unique_insert);
        uint64_t high_file_bytes = 0;
        for (const RootFileInfo& root_file : root_files)
        {
            high_file_bytes += root_file.file_size;
        }

        const uint64_t low_file_bytes = low_file_size_bytes<N>(options.tmp_dir, options.k_len);
        const uint64_t packed_low_kmer_bytes = packed_kmer_bytes_for_k<N>(options.k_len);
        ProgressPrinter progress(high_file_bytes + low_file_bytes);

        const uint64_t estimated_peak_memory =
            estimate_peak_memory_bytes<N>(expected_unique_insert, worker_count, hist_size);
        if (estimated_peak_memory > options.max_memory_bytes)
        {
            std::cerr << "insufficient-memory path is not implemented\n"
                << "estimated_peak_memory_bytes=" << estimated_peak_memory
                << " max_memory_bytes=" << options.max_memory_bytes << '\n';
            return 2;
        }

        AtomicHistogram global_histogram(hist_size);
        init_histogram(global_histogram);


        FlatConcurrentHashMap<N> hash_map(expected_unique_insert, worker_count);
        progress.start();

        {
            SPMCRingMemoryPool<HISTOGRAM_RING_CAPACITY> high_pool(HISTOGRAM_BLOCK_BYTES, 1);
            HighFrequencyInsertThreadPool<N, HISTOGRAM_RING_CAPACITY> high_threads(
                &high_pool,
                &hash_map,
                &global_histogram,
                worker_count,
                options.k_len,
                options.min_freq,
                options.max_freq,
                hist_size);

            high_threads.start();
            enqueue_high_records<N>(root_files, high_pool, &progress);
            high_threads.join();

            if (high_threads.insert_failed())
            {
                std::cerr << "failed to insert one or more high-frequency records" << std::endl;
                exit(-1);
            }
        }

        hash_map.seal();

        {
            SPMCRingMemoryPool<HISTOGRAM_RING_CAPACITY> low_pool(HISTOGRAM_BLOCK_BYTES, 1);
            LowFrequencyQueryThreadPool<N, HISTOGRAM_RING_CAPACITY> low_threads(
                &low_pool,
                &hash_map,
                &global_histogram,
                worker_count,
                options.k_len,
                options.min_freq,
                options.max_freq,
                hist_size);

            low_threads.start();
            enqueue_low_kmers<N>(options.k_len, low_pool, packed_low_kmer_bytes, &progress);
            low_threads.join();
        }

        write_histogram(options.output_file, global_histogram, options.min_freq);
        progress.finish();
        return 0;
    }

    int dispatch_by_k_len(const Options& options)
    {
        if (options.k_len <= 32)
        {
            return run_histogram_tool<1>(options);
        }
        else if (options.k_len <= 64)
        {
            return run_histogram_tool<2>(options);
        }
        else if (options.k_len <= MAX_K)
        {
            return run_histogram_tool<4>(options);
        }
        else
        {
            std::cerr << "unsupported k_len: " << options.k_len << '\n';
            exit(-1);
        }

    }
}

int main(int argc, char* argv[])
{
    const Options options = parse_options(argc, argv);
    return dispatch_by_k_len(options);

}
