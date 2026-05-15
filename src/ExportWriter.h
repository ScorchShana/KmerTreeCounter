#ifndef EXPORT_WRITER_HEADER
#define EXPORT_WRITER_HEADER

#include "definition.h"
#include "RingMemoryPool.h"
#include "kmer.h"

#include <cstdio>
#include <string>
#include <stdexcept>
#include <vector>
#include <array>
#include <thread>
#include <atomic>
#include <unordered_map>
#include <algorithm>

template <uint32_t N>
class ExportWriter
{

    constexpr static int MAX_SPIN_TIME = 1024;
    constexpr static int MAX_BACKOFF = 128;

private:
    // 连接到 FastqParser，用于获取含有低频 k-mer 的数据块
    RingMemoryPool<EXPORT_RING_MEMORY_POOL_CAPACITY> *ring_memory_pool = nullptr;
    // 为每个根节点（或者短前缀）维护一个文件句柄
    std::array<std::FILE *, EXPORT_FILES_SIZE> files{};
    // 每个根节点对应的最大缓冲容量（默认 8KB）
    static constexpr uint64_t ROOT_BUFFER_KMER_CAPACITY = EXPORT_ROOT_BUFFER_SIZE / sizeof(kmer<N>);
    // 对应每个根节点的本地缓存，当装满时通过 fwrite 顺序落盘
    std::array<std::vector<kmer<N>>, EXPORT_FILES_SIZE> root_buffers{};
    // 专门执行落盘操作的单线程
    std::thread worker_thread;
    std::atomic<bool> is_running{false};

public:
    ExportWriter() = delete;
    ExportWriter(const ExportWriter &) = delete;
    ExportWriter &operator=(const ExportWriter &) = delete;
    ExportWriter(ExportWriter &&) = delete;
    ExportWriter &operator=(ExportWriter &&) = delete;

    explicit ExportWriter(RingMemoryPool<EXPORT_RING_MEMORY_POOL_CAPACITY> *pool)
        : ring_memory_pool(pool)
    {
        for (uint64_t i = 0; i < files.size(); ++i)
        {
            files[i] = std::fopen((temp_dir + "low_" + std::to_string(i) + ".bin").c_str(), "wb");
            if (files[i] == nullptr) [[unlikely]]
            {
                throw std::runtime_error("Failed to open export file: " + temp_dir + "low_" + std::to_string(i) + ".bin");
            }
            root_buffers[i].reserve(ROOT_BUFFER_KMER_CAPACITY);
        }
    }

    ~ExportWriter()
    {
        for (uint64_t i = 0; i < files.size(); ++i)
        {
            if (files[i] != nullptr)
            {
                std::fclose(files[i]);
                files[i] = nullptr;
            }
        }
    }

    void start()
    {
        is_running = true;
        worker_thread = std::thread(&ExportWriter::writer_loop, this);
    }

    void join()
    {

        if (is_running && worker_thread.joinable())
        {
            worker_thread.join();
            is_running = false;
        }

        for (uint64_t i = 0; i < files.size(); ++i)
        {
            if (files[i] != nullptr)
            {
                std::fclose(files[i]);
                files[i] = nullptr;
            }
        }
    }

private:
    void writer_loop()
    {
        int spin_time = 0;
        int backoff = 1;

        char *block_ptr = nullptr;

        content_type content;

        while (true)
        {
            if (ring_memory_pool->consumer_try_dequeue(content))
            {
                block_ptr = content.data;

                ExportBlock<N> *export_kmer_block = reinterpret_cast<ExportBlock<N> *>(block_ptr);
                if (content.length > 0)
                {
                    process_block(export_kmer_block, content.length);
                }
                ring_memory_pool->consumer_enqueue(block_ptr);
            }
            else if (ring_memory_pool->producer_finished())
            {
                while (ring_memory_pool->consumer_try_dequeue(content))
                {
                    block_ptr = content.data;
                    ExportBlock<N> *export_kmer_block = reinterpret_cast<ExportBlock<N> *>(block_ptr);
                    if (content.length > 0)
                    {
                        process_block(export_kmer_block, content.length);
                    }
                    ring_memory_pool->consumer_enqueue(block_ptr);
                }
                break;
            }
            else
            {
                spin_time++;
                if (spin_time >= MAX_SPIN_TIME)
                {
                    std::this_thread::yield();
                    spin_time = 0;
                    backoff = 4;
                }
                else
                {
                    for (int i = 0; i < backoff; i++)
                    {
                        cpu_relax();
                    }
                    backoff = std::min(backoff * 2, MAX_BACKOFF);
                }
            }
        }

        flush_all_buffers();
    }

    void process_block(ExportBlock<N> *export_block_ptr, uint64_t kmer_amount)
    {
        for (uint64_t i = 0; i < kmer_amount; ++i)
        {
            const kmer<N> &val = export_block_ptr->k_mers[i];
            const uint64_t root_prefix = get_root_prefix(val);
            auto &buffer = root_buffers[root_prefix];
            buffer.push_back(val);
            if (buffer.size() >= ROOT_BUFFER_KMER_CAPACITY)
            {
                flush_root_buffer(root_prefix);
            }
        }
    }

    constexpr uint64_t get_root_prefix(const kmer<N> &val) const
    {
        constexpr uint32_t shift = 2 * (BASES_PER_U64T - ROOT_BASES);
        return val.data[0] >> shift;
    }

    void flush_root_buffer(const uint64_t root_prefix)
    {
        auto &buffer = root_buffers[root_prefix];

        const size_t data_size = buffer.size();
        if (std::fwrite(buffer.data(), sizeof(kmer<N>), data_size, files[root_prefix]) != data_size) [[unlikely]]
        {
            throw std::runtime_error("Failed to write k-mer data");
        }
        buffer.clear();
    }

    void flush_all_buffers()
    {
        for (uint64_t i = 0; i < root_buffers.size(); ++i)
        {
            if(!root_buffers[i].empty())
            {
                flush_root_buffer(i);
            }
        }
    }
};

#endif // EXPORT_WRITER_HEADER