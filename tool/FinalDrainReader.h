#ifndef FINAL_DRAIN_READER_HEADER
#define FINAL_DRAIN_READER_HEADER

#include "../src/definition.h"
#include "../src/kmer.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <aio.h>

template <uint32_t N>
class FinalDrainReader
{
    constexpr static uint64_t RAW_BUFFER_SIZE = 512 * 1024;
    constexpr static uint64_t RECORD_SIZE = sizeof(ExportRecord<N>);
    constexpr static uint64_t BUFFER_RECORD_CAPACITY = RAW_BUFFER_SIZE / RECORD_SIZE;
    constexpr static uint64_t BUFFER_BYTES = BUFFER_RECORD_CAPACITY * RECORD_SIZE;

    static_assert(BUFFER_RECORD_CAPACITY > 0, "ExportRecord is larger than FinalDrainReader raw buffer");

    int fd = -1;
    uint64_t record_amount = 0;
    uint64_t scheduled_record_count = 0;
    uint64_t delivered_record_count = 0;
    uint32_t current_buffer_index = 0;
    std::array<struct aiocb, 2> cbs{};
    std::array<char*, 2> buffer{};
    std::array<bool, 2> cbs_active{};
    std::array<uint64_t, 2> buffer_record_count{};
    std::array<uint64_t, 2> buffer_cursor{};

public:
    FinalDrainReader()
    {
        for (uint32_t i = 0; i < 2; ++i)
        {
            std::memset(&cbs[i], 0, sizeof(struct aiocb));
        }
    }

    FinalDrainReader(const FinalDrainReader&) = delete;
    FinalDrainReader& operator=(const FinalDrainReader&) = delete;
    FinalDrainReader(FinalDrainReader&&) = delete;
    FinalDrainReader& operator=(FinalDrainReader&&) = delete;

    ~FinalDrainReader()
    {
        close();
    }

    void open(const std::string& filename)
    {
        if (fd >= 0)
        {
            close();
        }

        fd = ::open(filename.c_str(), O_RDONLY);
        if (fd < 0) [[unlikely]]
        {
            std::cerr << "Failed to open file: " << filename << std::endl;
            std::exit(-1);
        }

        struct stat st;
        if (::fstat(fd, &st) != 0) [[unlikely]]
        {
            std::cerr << "Failed to stat file: " << filename << std::endl;
            std::exit(-1);
        }

        const uint64_t file_size = static_cast<uint64_t>(st.st_size);
        if (file_size % RECORD_SIZE != 0) [[unlikely]]
        {
            std::cerr << "Invalid final drain file size: " << file_size
                << " is not divisible by ExportRecord size "
                << RECORD_SIZE << std::endl;
            std::exit(-1);
        }

        ::posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);

        record_amount = file_size / RECORD_SIZE;
        scheduled_record_count = 0;
        delivered_record_count = 0;
        current_buffer_index = 0;

        for (uint32_t i = 0; i < 2; ++i)
        {
            std::memset(&cbs[i], 0, sizeof(struct aiocb));
            if (buffer[i] == nullptr)
            {
                buffer[i] = new char[BUFFER_BYTES];
            }
            cbs_active[i] = false;
            buffer_record_count[i] = 0;
            buffer_cursor[i] = 0;
        }

        async_read(0);
    }

    void close()
    {
        for (uint32_t i = 0; i < 2; ++i)
        {
            wait_read_finish(i);
        }

        if (fd >= 0)
        {
            ::close(fd);
            fd = -1;
        }

        for (uint32_t i = 0; i < 2; ++i)
        {
            delete[] buffer[i];
            buffer[i] = nullptr;
            cbs_active[i] = false;
            buffer_record_count[i] = 0;
            buffer_cursor[i] = 0;
            std::memset(&cbs[i], 0, sizeof(struct aiocb));
        }

        record_amount = 0;
        scheduled_record_count = 0;
        delivered_record_count = 0;
        current_buffer_index = 0;
    }

    uint64_t get_record_amount() const
    {
        return record_amount;
    }

    bool finished() const noexcept
    {
        return delivered_record_count >= record_amount;
    }

    uint64_t read_records(ExportRecord<N>* out_buffer, uint64_t max_records_to_read)
    {
        if (max_records_to_read == 0 || finished())
        {
            return 0;
        }

        uint64_t total_read = 0;
        while (total_read < max_records_to_read && !finished())
        {
            uint32_t idx = current_buffer_index;

            if (cbs_active[idx])
            {
                wait_read_finish(idx);
                async_read(1 - idx);
            }

            if (buffer_cursor[idx] >= buffer_record_count[idx])
            {
                current_buffer_index = 1 - idx;
                continue;
            }

            const uint64_t available = buffer_record_count[idx] - buffer_cursor[idx];
            if (available == 0)
            {
                current_buffer_index = 1 - idx;
                continue;
            }

            const uint64_t to_copy = std::min(max_records_to_read - total_read, available);
            const char* src = buffer[idx] + buffer_cursor[idx] * RECORD_SIZE;
            std::memcpy(out_buffer + total_read, src, static_cast<size_t>(to_copy * RECORD_SIZE));

            buffer_cursor[idx] += to_copy;
            delivered_record_count += to_copy;
            total_read += to_copy;

            if (buffer_cursor[idx] >= buffer_record_count[idx])
            {
                buffer_record_count[idx] = 0;
                buffer_cursor[idx] = 0;
                current_buffer_index = 1 - idx;
            }
        }

        return total_read;
    }

private:
    void async_read(uint32_t buffer_index)
    {
        if (scheduled_record_count >= record_amount || cbs_active[buffer_index]) [[unlikely]]
        {
            return;
        }

        const uint64_t remaining = record_amount - scheduled_record_count;
        const uint64_t read_count = std::min(BUFFER_RECORD_CAPACITY, remaining);
        if (read_count == 0)
        {
            return;
        }

        const uint64_t read_bytes = read_count * RECORD_SIZE;
        const uint64_t read_offset = scheduled_record_count * RECORD_SIZE;

        std::memset(&cbs[buffer_index], 0, sizeof(struct aiocb));
        cbs[buffer_index].aio_fildes = fd;
        cbs[buffer_index].aio_buf = buffer[buffer_index];
        cbs[buffer_index].aio_nbytes = static_cast<size_t>(read_bytes);
        cbs[buffer_index].aio_offset = static_cast<off_t>(read_offset);

        buffer_record_count[buffer_index] = read_count;
        buffer_cursor[buffer_index] = 0;

        if (::aio_read(&cbs[buffer_index]) != 0) [[unlikely]]
        {
            std::cerr << "aio_read failed" << std::endl;
            std::exit(-1);
        }

        scheduled_record_count += read_count;
        cbs_active[buffer_index] = true;
    }

    void wait_read_finish(uint32_t buffer_index)
    {
        if (!cbs_active[buffer_index])
        {
            return;
        }

        const aiocb* list[1] = { &cbs[buffer_index] };
        int err;
        while ((err = ::aio_error(&cbs[buffer_index])) == EINPROGRESS)
        {
            if (::aio_suspend(list, 1, nullptr) != 0 && errno != EINTR) [[unlikely]]
            {
                std::cerr << "aio_suspend failed" << std::endl;
                std::exit(-1);
            }
        }

        err = ::aio_error(&cbs[buffer_index]);
        if (err != 0) [[unlikely]]
        {
            std::cerr << "aio_read error" << std::endl;
            std::exit(-1);
        }

        const uint64_t expected_bytes = buffer_record_count[buffer_index] * RECORD_SIZE;
        const ssize_t n = ::aio_return(&cbs[buffer_index]);
        if (n != static_cast<ssize_t>(expected_bytes)) [[unlikely]]
        {
            std::cerr << "partial aio_read" << std::endl;
            std::exit(-1);
        }

        cbs_active[buffer_index] = false;
    }
};

#endif
