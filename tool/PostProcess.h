#ifndef POSTPROCESS_HEADER
#define POSTPROCESS_HEADER


#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <thread>
#include <algorithm>
#include <atomic>
#include <array>
#include <iostream>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <unordered_map>
#include <fstream>

inline const std::string temp_dir = "./project/KmerTreeCounter/build/tmp/";
inline const std::string output_dir = "./";
inline const uint64_t PAGE_SIZE = 4096;

template <uint32_t N>
using kmer = std::array<uint64_t, N>;

template <uint32_t N>
class PostProcess
{
    static constexpr uint64_t NUM_BUCKETS = 256;

    struct ExportRecord
    {
        kmer<N> key;
        uint32_t count;
    };

    std::vector<kmer<N>> highFreqKmer;
    std::vector<std::vector<ExportRecord>> buckets{NUM_BUCKETS};
    std::unordered_map<uint32_t, uint64_t> hist;   // count → num_kmers

    std::string output_low_filename;
    std::string output_root_filename;
    uint32_t num_threads;

public:
    PostProcess(const PostProcess &) = delete;
    PostProcess &operator=(const PostProcess &) = delete;
    PostProcess(PostProcess &&) = delete;
    PostProcess &operator=(PostProcess &&) = delete;

    explicit PostProcess(uint32_t threads)
        : num_threads(threads)
    {
        output_low_filename = output_dir + "final_low.bin";
        output_root_filename = output_dir + "final_root.bin";
        if (num_threads == 0)
            num_threads = 1;
    }

    void merge_all()
    {
        load_root_files();
        merge_buckets();
        process_low_bin();
        write_hist();
    }

    void write_hist()
    {
        if (hist.empty()) return;
        std::ofstream ofs(output_dir + "freq_hist.txt");
        if (!ofs) return;
        ofs << "count\tnum_kmers\n";
        std::vector<std::pair<uint32_t, uint64_t>> pairs(hist.begin(), hist.end());
        std::sort(pairs.begin(), pairs.end());
        for (const auto &p : pairs)
            ofs << p.first << '\t' << p.second << '\n';
    }

    void process_low_bin()
    {
        const std::string in_filename = temp_dir + "low.bin";
        const std::string tmp_prefix = temp_dir + ".low_sort_";

        { int fd = ::open(output_low_filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644); if (fd >= 0) ::close(fd); }

        struct stat st;
        if (::stat(in_filename.c_str(), &st) != 0) return;

        const uint64_t file_size = static_cast<uint64_t>(st.st_size);
        const uint64_t total_kmers = file_size / sizeof(kmer<N>);
        if (total_kmers == 0) return;

        uint64_t chunk_bytes = 256ULL * 1024 * 1024;
#ifdef _SC_AVPHYS_PAGES
        long pages = sysconf(_SC_AVPHYS_PAGES);
        long page_sz = sysconf(_SC_PAGESIZE);
        if (pages > 0 && page_sz > 0)
        {
            uint64_t avail_mem = static_cast<uint64_t>(pages) * static_cast<uint64_t>(page_sz);
            chunk_bytes = avail_mem / 3 / (static_cast<uint64_t>(num_threads) + 1);
            chunk_bytes = std::max<uint64_t>(chunk_bytes, 64ULL << 20);
            chunk_bytes = std::min<uint64_t>(chunk_bytes, 1024ULL << 20);
            chunk_bytes &= ~(PAGE_SIZE - 1);
        }
#endif
        const uint64_t CHUNK_KMERS = chunk_bytes / sizeof(kmer<N>);
        const uint64_t num_chunks = (total_kmers + CHUNK_KMERS - 1) / CHUNK_KMERS;

        std::vector<std::string> chunk_files(num_chunks);
        std::atomic<uint64_t> next_chunk{0};
        std::vector<std::thread> sort_threads;
        sort_threads.reserve(std::min(num_chunks, static_cast<uint64_t>(num_threads)));

        for (uint32_t t = 0; t < num_threads && t < num_chunks; ++t)
        {
            sort_threads.emplace_back([&]()
            {
                while (true)
                {
                    uint64_t chunk_id = next_chunk.fetch_add(1, std::memory_order_relaxed);
                    if (chunk_id >= num_chunks) break;
                    sort_low_chunk(in_filename, chunk_id, CHUNK_KMERS, total_kmers, tmp_prefix, chunk_files);
                }
            });
        }
        for (auto &th : sort_threads) if (th.joinable()) th.join();

        merge_sorted_chunks(chunk_files, output_low_filename);

        for (auto &f : chunk_files) if (!f.empty()) ::unlink(f.c_str());

        // 低频 k-mer 加入直方图 
        if (::stat(output_low_filename.c_str(), &st) == 0 && st.st_size > 0)
            hist[1] += static_cast<uint64_t>(st.st_size) / sizeof(kmer<N>);
    }

private:
    void load_root_files()
    {
        for (uint64_t i = 0; i < NUM_BUCKETS; ++i)
        {
            const std::string filename = temp_dir + "root_" + std::to_string(i) + ".bin";
            struct stat st;
            if (::stat(filename.c_str(), &st) != 0) continue;

            const uint64_t file_size = static_cast<uint64_t>(st.st_size);
            if (file_size == 0) continue;

            const uint64_t record_count = file_size / sizeof(ExportRecord);
            buckets[i].reserve(buckets[i].size() + record_count);
            highFreqKmer.reserve(highFreqKmer.size() + record_count);

            int fd = ::open(filename.c_str(), O_RDONLY);
            if (fd < 0) [[unlikely]] continue;

            std::vector<ExportRecord> buf(record_count);
            const ssize_t bytes = ::read(fd, buf.data(), file_size);
            if (bytes != static_cast<ssize_t>(file_size)) [[unlikely]]
                std::cerr << "[PostProcess] Short read on " << filename << std::endl;
            ::close(fd);

            buckets[i].insert(buckets[i].end(), buf.begin(), buf.begin() + record_count);
            for (auto &rec : buf) {
                highFreqKmer.push_back(rec.key);
                // Bloom Filter 首次误判为低频 → 计数少了 1, 补偿回去
                hist[static_cast<uint32_t>(rec.count + 1)]++;
            }
        }
        std::sort(highFreqKmer.begin(), highFreqKmer.end());
        highFreqKmer.erase(std::unique(highFreqKmer.begin(), highFreqKmer.end()), highFreqKmer.end());
    }

    void merge_buckets()
    {
        int out_fd = ::open(output_root_filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (out_fd < 0) [[unlikely]] { std::cerr << "[PostProcess] Failed to open output\n"; std::exit(-1); }

        constexpr uint64_t WRITE_BUF_SIZE = 2ULL * 1024 * 1024;
        std::vector<char> write_buf(WRITE_BUF_SIZE);
        uint64_t buf_offset = 0;

        for (uint64_t i = 0; i < NUM_BUCKETS; ++i)
        {
            const auto &bucket = buckets[i];
            if (bucket.empty()) continue;

            const char *data = reinterpret_cast<const char *>(bucket.data());
            uint64_t remaining = bucket.size() * sizeof(ExportRecord);
            uint64_t offset = 0;

            while (remaining > 0)
            {
                const uint64_t space = WRITE_BUF_SIZE - buf_offset;
                const uint64_t to_copy = std::min(space, remaining);
                std::memcpy(write_buf.data() + buf_offset, data + offset, to_copy);
                buf_offset += to_copy;
                offset += to_copy;
                remaining -= to_copy;
                if (buf_offset >= WRITE_BUF_SIZE) {
                    ::write(out_fd, write_buf.data(), WRITE_BUF_SIZE);
                    buf_offset = 0;
                }
            }
        }
        if (buf_offset > 0) ::write(out_fd, write_buf.data(), buf_offset);
        ::fsync(out_fd); ::close(out_fd);
    }



    //外部排序

    void sort_low_chunk(const std::string &in_filename, uint64_t chunk_id,
                        uint64_t chunk_capacity, uint64_t total_kmers,
                        const std::string &tmp_prefix, std::vector<std::string> &out_files)
    {
        uint64_t offset = chunk_id * chunk_capacity;
        uint64_t count = std::min(chunk_capacity, total_kmers - offset);
        std::vector<kmer<N>> buf(count);

        int fd = ::open(in_filename.c_str(), O_RDONLY);
        if (fd < 0) [[unlikely]] { return; }
        ssize_t bytes = ::pread(fd, buf.data(), count * sizeof(kmer<N>), static_cast<off_t>(offset * sizeof(kmer<N>)));
        ::close(fd);

        if (bytes != static_cast<ssize_t>(count * sizeof(kmer<N>))) [[unlikely]] {
            count = (bytes > 0) ? static_cast<uint64_t>(bytes) / sizeof(kmer<N>) : 0;
            if (count == 0) return;
            buf.resize(count);
        }
        std::sort(buf.begin(), buf.end());

        std::string tmp_name = tmp_prefix + std::to_string(chunk_id) + ".bin";
        int out_fd = ::open(tmp_name.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (out_fd < 0) [[unlikely]] { return; }
        ::write(out_fd, buf.data(), count * sizeof(kmer<N>));
        ::fsync(out_fd); ::close(out_fd);
        out_files[chunk_id] = tmp_name;
    }

    void merge_sorted_chunks(const std::vector<std::string> &chunk_files, const std::string &out_filename)
    {
        struct ChunkInfo { std::string fname; uint64_t total; };
        std::vector<ChunkInfo> valid_chunks;
        for (auto &f : chunk_files) {
            if (f.empty()) continue;
            struct stat st; if (::stat(f.c_str(), &st) != 0) continue;
            valid_chunks.push_back({f, static_cast<uint64_t>(st.st_size) / sizeof(kmer<N>)});
        }
        if (valid_chunks.empty()) return;

        const uint64_t num_ways = valid_chunks.size();
        std::vector<int> fds(num_ways, -1);
        for (uint64_t i = 0; i < num_ways; ++i) {
            fds[i] = ::open(valid_chunks[i].fname.c_str(), O_RDONLY);
            if (fds[i] < 0) [[unlikely]] { valid_chunks[i].total = 0; }
        }

        static constexpr uint64_t READAHEAD_KMERS = 8192;
        std::vector<std::vector<kmer<N>>> buffers(num_ways);
        std::vector<uint64_t> buf_pos(num_ways, 0), buf_end(num_ways, 0), remaining(num_ways);
        for (uint64_t i = 0; i < num_ways; ++i) {
            remaining[i] = valid_chunks[i].total;
            buffers[i].resize(READAHEAD_KMERS);
            refill_buffer(i, fds[i], buffers[i], buf_pos[i], buf_end[i], remaining[i]);
        }

        constexpr uint64_t WRITE_BUF_BYTES = 2ULL * 1024 * 1024;
        int out_fd = ::open(out_filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (out_fd < 0) [[unlikely]] { for (auto fd : fds) if (fd >= 0) ::close(fd); return; }
        std::vector<char> out_buf(WRITE_BUF_BYTES);
        uint64_t out_offset = 0;

        struct HeapItem { kmer<N> key; uint64_t way; bool operator<(const HeapItem &o) const { return o.key < key; } };
        std::vector<HeapItem> heap; heap.reserve(num_ways);
        for (uint64_t i = 0; i < num_ways; ++i)
            if (buf_end[i] > buf_pos[i] || remaining[i] > 0)
                heap.push_back({current_kmer(i, buffers[i], buf_pos[i], buf_end[i], fds[i], remaining[i]), i});
        std::make_heap(heap.begin(), heap.end());

        kmer<N> *last_written = nullptr;

        while (!heap.empty())
        {
            std::pop_heap(heap.begin(), heap.end());
            HeapItem top = heap.back(); heap.pop_back();

            if (std::binary_search(highFreqKmer.begin(), highFreqKmer.end(), top.key)) {
                uint64_t w = top.way;
                buf_pos[w]++;
                if (buf_pos[w] >= buf_end[w] && remaining[w] > 0)
                    refill_buffer(w, fds[w], buffers[w], buf_pos[w], buf_end[w], remaining[w]);
                if (buf_pos[w] < buf_end[w] || remaining[w] > 0) {
                    heap.push_back({current_kmer(w, buffers[w], buf_pos[w], buf_end[w], fds[w], remaining[w]), w});
                    std::push_heap(heap.begin(), heap.end());
                }
                continue;
            }

            if (last_written == nullptr || !(top.key == *last_written)) {
                constexpr uint64_t ksz = sizeof(kmer<N>);
                if (out_offset + ksz > WRITE_BUF_BYTES) { ::write(out_fd, out_buf.data(), out_offset); out_offset = 0; }
                std::memcpy(out_buf.data() + out_offset, &top.key, ksz);
                last_written = reinterpret_cast<kmer<N>*>(out_buf.data() + out_offset);
                out_offset += ksz;
            }

            uint64_t w = top.way; buf_pos[w]++;
            if (buf_pos[w] >= buf_end[w] && remaining[w] > 0)
                refill_buffer(w, fds[w], buffers[w], buf_pos[w], buf_end[w], remaining[w]);
            if (buf_pos[w] < buf_end[w] || remaining[w] > 0) {
                heap.push_back({current_kmer(w, buffers[w], buf_pos[w], buf_end[w], fds[w], remaining[w]), w});
                std::push_heap(heap.begin(), heap.end());
            }
        }

        if (out_offset > 0) ::write(out_fd, out_buf.data(), out_offset);
        ::fsync(out_fd); ::close(out_fd);
        for (auto fd : fds) if (fd >= 0) ::close(fd);
    }

    static kmer<N>& current_kmer(uint64_t way, std::vector<kmer<N>> &buf,
                                  uint64_t &pos, uint64_t &end, int fd, uint64_t &remaining)
    {
        if (pos >= end) refill_buffer(way, fd, buf, pos, end, remaining);
        return buf[pos];
    }

    static void refill_buffer(uint64_t way, int fd, std::vector<kmer<N>> &buf,
                              uint64_t &pos, uint64_t &end, uint64_t &remaining)
    {
        if (remaining == 0 || fd < 0) { pos = end = 0; return; }
        uint64_t to_read = std::min(static_cast<uint64_t>(buf.size()), remaining);
        ssize_t bytes = ::read(fd, buf.data(), to_read * sizeof(kmer<N>));
        if (bytes <= 0) { pos = end = 0; remaining = 0; return; }
        pos = 0;
        end = static_cast<uint64_t>(bytes) / sizeof(kmer<N>);
        remaining -= end;
    }
};

#endif