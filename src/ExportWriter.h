#ifndef EXPORT_WRITER_HEADER
#define EXPORT_WRITER_HEADER

#include "definition.h"
#include "RingMemoryPool.h"
#include "kmer.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <array>
#include <thread>
#include <atomic>
#include <unordered_map>
#include <algorithm>

#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <errno.h>


template <uint32_t N>
class ExportWriter
{

    constexpr static int MAX_SPIN_TIME = 1024;
    constexpr static int MAX_BACKOFF = 128;

    constexpr static uint64_t RAW_BUFFER_SIZE = 4ULL * 1024 * 1024; // 4MB 原始数据块大小
    constexpr static uint64_t BUFFER_KMER_CAPACITY = RAW_BUFFER_SIZE / sizeof(kmer<N>);

private:
    // 连接到 FastqParser，用于获取含有低频 k-mer 的数据块
    RingMemoryPool<EXPORT_RING_MEMORY_POOL_CAPACITY>* ring_memory_pool = nullptr;

    int fd;
    kmer<N>* buffer;
    uint64_t buffer_count;

    // 专门执行落盘操作的单线程
    std::thread worker_thread;
    std::atomic<bool> is_running{ false };

public:
    ExportWriter() = delete;
    ExportWriter(const ExportWriter&) = delete;
    ExportWriter& operator=(const ExportWriter&) = delete;
    ExportWriter(ExportWriter&&) = delete;
    ExportWriter& operator=(ExportWriter&&) = delete;

    explicit ExportWriter(RingMemoryPool<EXPORT_RING_MEMORY_POOL_CAPACITY>* pool)
        : ring_memory_pool(pool), fd(-1), buffer(nullptr), buffer_count(0)
    {
        fd = ::open((temp_dir + "low.bin").c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (fd < 0) [[unlikely]]
        {
            std::cerr << "Failed to open export file: " << temp_dir + "low.bin" << std::endl;
            std::exit(-1);
        }
    }

    void join()
    {
        if (is_running && worker_thread.joinable())
        {
            worker_thread.join();
            is_running = false;
        }

        const std::string filename = temp_dir + "low.bin";

        if (fd >= 0)
        {
            ::close(fd);
            fd = -1;
        }
    }

    ~ExportWriter()
    {
        if (fd >= 0)
        {
            ::close(fd);
            fd = -1;
        }
        if (buffer != nullptr)
        {
            delete[] buffer;
            buffer = nullptr;
        }
    }

    void start()
    {
        is_running = true;
        worker_thread = std::thread(&ExportWriter::writer_loop, this);
    }

private:
    void writer_loop()
    {

        buffer = new kmer<N>[BUFFER_KMER_CAPACITY];

        int spin_time = 0;
        int backoff = 1;

        char* block_ptr = nullptr;

        content_type content;

        while (true)
        {
            if (ring_memory_pool->consumer_try_dequeue(content))
            {
                block_ptr = content.data;

                ExportBlock<N>* export_kmer_block = reinterpret_cast<ExportBlock<N> *>(block_ptr);
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
                    ExportBlock<N>* export_kmer_block = reinterpret_cast<ExportBlock<N> *>(block_ptr);
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

        flush_buffer();

        delete[] buffer;
        buffer = nullptr;
    }

    void process_block(ExportBlock<N>* export_block_ptr, uint64_t kmer_amount)
    {
        if (buffer_count + kmer_amount >= BUFFER_KMER_CAPACITY) [[unlikely]]
        {
            uint64_t to_write = BUFFER_KMER_CAPACITY - buffer_count;
            std::memcpy(buffer + buffer_count, export_block_ptr->k_mers.data(), to_write * sizeof(kmer<N>));
            flush_buffer();
            uint64_t remaining = kmer_amount - to_write;
            std::memcpy(buffer, export_block_ptr->k_mers.data() + to_write, remaining * sizeof(kmer<N>));
            buffer_count = remaining;
        }
        else
        {
            std::memcpy(buffer + buffer_count, export_block_ptr->k_mers.data(), kmer_amount * sizeof(kmer<N>));
            buffer_count += kmer_amount;
        }
    }

    constexpr uint64_t get_root_prefix(const kmer<N>& val) const
    {
        constexpr uint32_t shift = 2 * (BASES_PER_U64T - ROOT_BASES);
        return val.data[0] >> shift;
    }

    void flush_buffer()
    {
        if (!write_all(buffer, sizeof(kmer<N>) * buffer_count)) [[unlikely]]
        {
            std::cerr << "Failed to write k-mer data" << std::endl;
            std::exit(-1);
        }
        buffer_count = 0;
    }

    bool write_all(void* data, size_t count) {
        const uint8_t* p = static_cast<const uint8_t*>(data);
        size_t left = count;

        while (left > 0) {
            ssize_t n = ::write(fd, p, left);

            if (n > 0) [[unlikely]]
            {
                p += n;
                left -= static_cast<size_t>(n);
                continue;
            }
            else if (n < 0 && errno == EINTR) [[unlikely]]
            {
                continue;
            }

            return false;
        }

        return true;
    }
};

#endif // EXPORT_WRITER_HEADER
