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
#include <aio.h>


template <uint32_t N>
class ExportWriter
{

    constexpr static int MAX_SPIN_TIME = 1024;
    constexpr static int MAX_BACKOFF = 128;

    constexpr static uint64_t RAW_BUFFER_SIZE = 1ULL * 1024 * 1024; // 1MB 原始数据块大小
    constexpr static uint64_t BUFFER_KMER_CAPACITY = RAW_BUFFER_SIZE / sizeof(kmer<N>);

private:
    // 连接到 FastqParser，用于获取含有低频 k-mer 的数据块
    RingMemoryPool<EXPORT_RING_MEMORY_POOL_CAPACITY>* ring_memory_pool = nullptr;

    struct aiocb cbs[2];
    int fd;
    bool cbs_active[2];
    std::array<kmer<N>*, 2> buffer{};
    std::array<uint32_t, 2> buffer_count{ 0, 0 }; // 当前缓冲区中 k-mer 的数量
    uint32_t current_buffer_index = 0;
    uint64_t file_offset = 0; // 已经写入文件的字节数

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
        : ring_memory_pool(pool), fd(-1)
    {
        fd = ::open((temp_dir + "low.bin").c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (fd < 0) [[unlikely]]
        {
            std::cerr << "Failed to open export file: " << temp_dir + "low.bin" << std::endl;
            std::exit(-1);
        }

        for (int i = 0; i < 2; ++i)
        {
            std::memset(&cbs[i], 0, sizeof(struct aiocb));
            buffer[i] = nullptr;
            buffer_count[i] = 0;
            cbs_active[i] = false;
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
        for (auto& buf : buffer) {
            if (buf != nullptr)
            {
                delete[] buf;
                buf = nullptr;
            }
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

        std::memset(cbs, 0, sizeof(cbs));
        buffer[0] = new kmer<N>[BUFFER_KMER_CAPACITY];
        buffer[1] = new kmer<N>[BUFFER_KMER_CAPACITY];

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
                if (content.length > 0) [[likely]]
                {
                    process_block(export_kmer_block, content.length);
                }
                else {
                    ring_memory_pool->consumer_enqueue(block_ptr);
                }

            }
            else if (ring_memory_pool->producer_finished())
            {
                while (ring_memory_pool->consumer_try_dequeue(content))
                {
                    block_ptr = content.data;
                    ExportBlock<N>* export_kmer_block = reinterpret_cast<ExportBlock<N> *>(block_ptr);
                    if (content.length > 0) [[likely]]
                    {
                        process_block(export_kmer_block, content.length);
                    }
                    else {
                        ring_memory_pool->consumer_enqueue(block_ptr);
                    }

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


        wait_write_finish(1 - current_buffer_index);

        if (buffer_count[current_buffer_index] > 0)
        {
            async_write(current_buffer_index);
            wait_write_finish(current_buffer_index);
        }

        delete[] buffer[0];
        delete[] buffer[1];
        buffer[0] = nullptr;
        buffer[1] = nullptr;
    }

    void process_block(ExportBlock<N>* export_block_ptr, uint64_t kmer_amount)
    {
        if (buffer_count[current_buffer_index] + kmer_amount >= BUFFER_KMER_CAPACITY) [[unlikely]]
        {
            uint64_t to_write = BUFFER_KMER_CAPACITY - buffer_count[current_buffer_index];
            std::memcpy(buffer[current_buffer_index] + buffer_count[current_buffer_index], export_block_ptr->k_mers.data(), to_write * sizeof(kmer<N>));
            uint64_t remaining = kmer_amount - to_write;

            async_write(current_buffer_index);

            current_buffer_index = 1 - current_buffer_index;
            wait_write_finish(current_buffer_index);


            std::memcpy(buffer[current_buffer_index] + buffer_count[current_buffer_index], export_block_ptr->k_mers.data() + to_write, remaining * sizeof(kmer<N>));
            buffer_count[current_buffer_index] = remaining;

            ring_memory_pool->consumer_enqueue(reinterpret_cast<char*>(export_block_ptr));

        }
        else
        {
            std::memcpy(buffer[current_buffer_index] + buffer_count[current_buffer_index], export_block_ptr->k_mers.data(), kmer_amount * sizeof(kmer<N>));
            buffer_count[current_buffer_index] += kmer_amount;
            ring_memory_pool->consumer_enqueue(reinterpret_cast<char*>(export_block_ptr));
        }
    }

    void wait_write_finish(const uint32_t buffer_index)
    {
        if (!cbs_active[buffer_index]) [[unlikely]]
        {
            return;
        }

        const aiocb* list[1] = { &cbs[buffer_index] };
        int err;
        while ((err = ::aio_error(&cbs[buffer_index])) == EINPROGRESS)
        {
            if (::aio_suspend(list, 1, nullptr) != 0 && errno != EINTR) [[unlikely]]
            {
                std::cerr << "aio_suspend failed\n";
                std::exit(-1);
            }
        }

        err = ::aio_error(&cbs[buffer_index]);
        if (err != 0) [[unlikely]]
        {
            std::cerr << "aio_write error\n";
            std::exit(-1);
        }

        ssize_t n = ::aio_return(&cbs[buffer_index]);
        if (n != static_cast<ssize_t>(cbs[buffer_index].aio_nbytes)) [[unlikely]]
        {
            std::cerr << "partial aio_write\n";
            std::exit(-1);
        }

        cbs_active[buffer_index] = false;
    }

    void async_write(const uint32_t buffer_index)
    {
        std::memset(&cbs[buffer_index], 0, sizeof(struct aiocb));
        cbs[buffer_index].aio_fildes = fd;
        cbs[buffer_index].aio_buf = buffer[buffer_index];
        cbs[buffer_index].aio_nbytes = buffer_count[buffer_index] * sizeof(kmer<N>);
        cbs[buffer_index].aio_offset = file_offset;

        if (::aio_write(&cbs[buffer_index]) != 0) [[unlikely]]
        {
            std::cerr << "aio_write failed\n";
            std::exit(-1);
        }

        file_offset += buffer_count[buffer_index] * sizeof(kmer<N>);
        buffer_count[buffer_index] = 0;
        cbs_active[buffer_index] = true;
    }

    constexpr uint64_t get_root_prefix(const kmer<N>& val) const
    {
        constexpr uint32_t shift = 2 * (BASES_PER_U64T - ROOT_BASES);
        return val.data[0] >> shift;
    }

};

#endif // EXPORT_WRITER_HEADER
