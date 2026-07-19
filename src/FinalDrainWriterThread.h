#ifndef FINAL_DRAIN_WRITER_THREAD_HEADER
#define FINAL_DRAIN_WRITER_THREAD_HEADER

#include "definition.h"
#include "RingMemoryPool.h"
#include "SpinBackoff.h"

#include <aio.h>
#include <fcntl.h>
#include <unistd.h>

#include <thread>
#include <array>
#include <cstring>
#include <iostream>

class FinalDrainWriterThread
{
    static constexpr uint64_t AIO_BUFFER_SIZE = 512 * 1024;

    RingMemoryPool<FINAL_DRAIN_RING_POOL_CAPACITY> pool_;
    int fd_;
    std::thread thread_;

    struct aiocb cbs_[2];
    bool cbs_active_[2];
    std::array<char*, 2> buffer_;
    std::array<uint64_t, 2> buffer_count_;
    uint32_t current_buffer_index_ = 0;
    uint64_t file_offset_ = 0;

public:
    FinalDrainWriterThread(uint32_t block_size, uint32_t producer_count)
        : pool_(block_size, producer_count), fd_(-1) {}

    auto* pool() { return &pool_; }

    void start()
    {
        fd_ = ::open((temp_dir + "high.bin").c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (fd_ < 0) [[unlikely]]
        {
            std::cerr << "Failed to open high.bin\n";
            std::exit(-1);
        }
        cbs_active_[0] = cbs_active_[1] = false;
        buffer_[0] = buffer_[1] = nullptr;
        buffer_count_[0] = buffer_count_[1] = 0;
        thread_ = std::thread(&FinalDrainWriterThread::writer_loop, this);
    }

    void join()
    {
        if (thread_.joinable()) thread_.join();
        if (fd_ >= 0)
        {
            ::close(fd_);
            fd_ = -1;
        }
        for (int i = 0; i < 2; i++)
        {
            if (buffer_[i] != nullptr)
            {
                delete[] buffer_[i];
                buffer_[i] = nullptr;
            }
        }
    }

    FinalDrainWriterThread(const FinalDrainWriterThread&) = delete;
    FinalDrainWriterThread& operator=(const FinalDrainWriterThread&) = delete;

    ~FinalDrainWriterThread()
    {
        join();
    }

private:
    void writer_loop()
    {
        buffer_[0] = new char[AIO_BUFFER_SIZE]();
        buffer_[1] = new char[AIO_BUFFER_SIZE]();

        SpinBackoff<128, 128, 256 * 1024> backoff;
        content_type content;

        while (true)
        {
            if (pool_.consumer_try_dequeue(content))
            {
                backoff.decay();
                process_block(content);
            }
            else if (pool_.producer_finished())
            {
                while (pool_.consumer_try_dequeue(content))
                {
                    process_block(content);
                    cpu_relax();
                }
                break;
            }
            else
            {
                backoff.backoff();
            }
        }

        wait_write_finish(1 - current_buffer_index_);
        if (buffer_count_[current_buffer_index_] > 0)
        {
            async_write(current_buffer_index_);
            wait_write_finish(current_buffer_index_);
        }

        delete[] buffer_[0]; buffer_[0] = nullptr;
        delete[] buffer_[1]; buffer_[1] = nullptr;
    }

    void process_block(const content_type& content)
    {
        if (content.length == 0) [[unlikely]]
        {
            pool_.consumer_enqueue(content.data);
            return;
        }

        if (buffer_count_[current_buffer_index_] + content.length >= AIO_BUFFER_SIZE)
        {
            uint64_t to_copy = AIO_BUFFER_SIZE - buffer_count_[current_buffer_index_];
            std::memcpy(buffer_[current_buffer_index_] + buffer_count_[current_buffer_index_],
                content.data, to_copy);
            uint64_t remaining = content.length - to_copy;
            buffer_count_[current_buffer_index_] += to_copy;

            async_write(current_buffer_index_);
            current_buffer_index_ = 1 - current_buffer_index_;
            wait_write_finish(current_buffer_index_);

            std::memcpy(buffer_[current_buffer_index_], content.data + to_copy, remaining);
            buffer_count_[current_buffer_index_] = remaining;
        }
        else
        {
            std::memcpy(buffer_[current_buffer_index_] + buffer_count_[current_buffer_index_],
                content.data, content.length);
            buffer_count_[current_buffer_index_] += content.length;
        }
        pool_.consumer_enqueue(content.data);
    }

    void async_write(uint32_t idx)
    {
        std::memset(&cbs_[idx], 0, sizeof(struct aiocb));
        cbs_[idx].aio_fildes = fd_;
        cbs_[idx].aio_buf = buffer_[idx];
        cbs_[idx].aio_nbytes = buffer_count_[idx];
        cbs_[idx].aio_offset = file_offset_;

        if (::aio_write(&cbs_[idx]) != 0) [[unlikely]]
        {
            std::cerr << "aio_write failed\n";
            std::exit(-1);
        }

        file_offset_ += buffer_count_[idx];
        buffer_count_[idx] = 0;
        cbs_active_[idx] = true;
    }

    void wait_write_finish(uint32_t idx)
    {
        if (!cbs_active_[idx]) [[unlikely]]
        {
            return;
        }

        const aiocb* list[1] = { &cbs_[idx] };
        int err;
        while ((err = ::aio_error(&cbs_[idx])) == EINPROGRESS)
        {
            if (::aio_suspend(list, 1, nullptr) != 0 && errno != EINTR) [[unlikely]]
            {
                std::cerr << "aio_suspend failed\n";
                std::exit(-1);
            }
        }

        err = ::aio_error(&cbs_[idx]);
        if (err != 0) [[unlikely]]
        {
            std::cerr << "aio_write error\n";
            std::exit(-1);
        }

        ssize_t n = ::aio_return(&cbs_[idx]);
        if (n != static_cast<ssize_t>(cbs_[idx].aio_nbytes)) [[unlikely]]
        {
            std::cerr << "partial aio_write\n";
            std::exit(-1);
        }

        cbs_active_[idx] = false;
    }
};

#endif // FINAL_DRAIN_WRITER_THREAD_HEADER
