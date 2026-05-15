#ifndef FASTQ_READER_HEADER
#define FASTQ_READER_HEADER

#include "definition.h"
#include "RingMemoryPool.h"

#include <cassert>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <utility>
#include <iostream>

template <uint32_t N>
class FastqReader
{

    using content_type = std::pair<char *, uint64_t>;

public:
    int k;

    explicit FastqReader(const std::string &filename, const int in_k, uint64_t chunk_size, RingMemoryPool<RING_MEMORY_POOL_CAPACITY> *in_ring_memory_pool_ptr)
        : filename_(filename), ring_memory_pool_ptr_(in_ring_memory_pool_ptr), file_buffer_(chunk_size),
          k(in_k), chunk_size_(chunk_size)
    {
        assert(ring_memory_pool_ptr_ != nullptr && "RingMemoryPool pointer must not be null");
        assert(k > 0 && "k must be positive");
        assert(k < 128 && "k must be < 128");
        assert(ring_memory_pool_ptr_->blockSize() > 1024 && "block_size must be > 1024");

        if (ring_memory_pool_ptr_ == nullptr)
        {
            throw std::invalid_argument("RingMemoryPool pointer is null");
        }
        fd_ = ::open(filename.data(), O_RDONLY | O_LARGEFILE);
        if (fd_ == -1)
        {
            throw std::runtime_error("Failed to open file");
        }

        posix_fadvise(fd_, 0, 0, POSIX_FADV_SEQUENTIAL); // 顺序访问提示

        struct stat st;
        if (fstat(fd_, &st) == -1)
        {
            throw std::runtime_error("Failed to get file size");
        }

        file_size_ = st.st_size;
    }

    ~FastqReader()
    {
        if (fd_ != -1)
        {
            ::close(fd_);
        }
    }

    void read()
    {
        assert(ring_memory_pool_ptr_ != nullptr && "RingMemoryPool pointer must not be null");
        char *block_ptr = nullptr;
        const uint64_t block_size = ring_memory_pool_ptr_->blockSize();
        assert(k > 0 && "k must be positive");
        assert(k < 128 && "k must be < 128");
        assert(block_size > 1024 && "block_size must be > 1024");

        const uint64_t overlap = (k > 1) ? static_cast<uint64_t>(k - 1) : 0;
        uint64_t write_size = 0;
        uint64_t last_newline_pos = kNoNewlineInBlock;
        bool has_block = false;

        uint64_t input_pos = 0;
        uint64_t input_size = 0;
        bool eof = false;
        bool stop = false;

#ifdef TEST_MODE
        double finish_percent = 0.0;
        double last_reported_percent = 0.0;
#endif

        while (!stop)
        {
            if (!has_block)
            {
                if (!acquire_block(block_ptr, write_size, last_newline_pos, has_block))
                {
                    break;
                }
            }

            if (input_pos >= input_size && !eof)
            {
                const ssize_t bytes_read = ::read(fd_, file_buffer_.data(), chunk_size_);
                if (bytes_read < 0) [[unlikely]]
                {
                    throw std::runtime_error("Failed to read fastq data");
                }

                if (bytes_read == 0)
                {
                    eof = true;
                }
                else
                {
                    input_pos = 0;
                    input_size = static_cast<uint64_t>(bytes_read);
                    have_read_ += bytes_read;
#ifdef TEST_MODE
                    last_reported_percent = finish_percent;
                    finish_percent = static_cast<double>(have_read_) / static_cast<double>(file_size_);
                    if (finish_percent - last_reported_percent >= 0.01 || finish_percent == 1.0)
                    {
                        std::cout << "FastqReader progress: " << static_cast<int>(finish_percent * 100) << "%\n";
                    }
#endif

                    if (have_read_ + static_cast<off_t>(chunk_size_) < file_size_)
                    {
                        ::posix_fadvise(fd_, have_read_, static_cast<size_t>(chunk_size_), POSIX_FADV_WILLNEED);
                    }
                }
            }

            if (input_pos >= input_size)
            {
                if (eof)
                {
                    publish_current_block(block_ptr, write_size, last_newline_pos, has_block);
                    break;
                }
                continue;
            }

            const char *input_begin = file_buffer_.data();

            if (state_ != State::ReadSequence)
            {
                const char *cur = input_begin + input_pos;
                const uint64_t remain = input_size - input_pos;
                const void *newline_ptr = std::memchr(cur, '\n', remain);

                if (newline_ptr == nullptr)
                {
                    input_pos = input_size;
                    continue;
                }

                const char *newline = static_cast<const char *>(newline_ptr);
                input_pos = static_cast<uint64_t>(newline - input_begin) + 1;
                state_ = advance_state_on_newline(state_);
                continue;
            }

            // HOT PATH - sequence line is copied in batches instead of char-by-char.
            const char *seq_cur = input_begin + input_pos;
            const uint64_t seq_remain = input_size - input_pos;
            const void *newline_ptr = std::memchr(seq_cur, '\n', seq_remain);
            const uint64_t sequence_len = (newline_ptr == nullptr)
                                              ? seq_remain
                                              : static_cast<uint64_t>(static_cast<const char *>(newline_ptr) - seq_cur);

            uint64_t copied = 0;
            while (copied < sequence_len)
            {
                if (write_size == block_size)
                {
                    store_overlap_from_block_end(block_ptr, write_size, last_newline_pos, overlap);
                    publish_current_block(block_ptr, write_size, last_newline_pos, has_block);
                    if (!acquire_block(block_ptr, write_size, last_newline_pos, has_block))
                    {
                        stop = true;
                        break;
                    }
                }

                const uint64_t block_remain = block_size - write_size;
                const uint64_t to_copy = (sequence_len - copied < block_remain) ? (sequence_len - copied) : block_remain;
                std::memcpy(block_ptr + write_size, seq_cur + copied, to_copy);
                write_size += to_copy;
                copied += to_copy;
                input_pos += to_copy;
            }

            if (stop)
            {
                break;
            }

            if (newline_ptr == nullptr)
            {
                continue;
            }

            if (write_size == block_size)
            {
                left_buffer_size_ = 0;
                publish_current_block(block_ptr, write_size, last_newline_pos, has_block);
                if (!acquire_block(block_ptr, write_size, last_newline_pos, has_block))
                {
                    break;
                }
            }

            last_newline_pos = write_size;
            block_ptr[write_size++] = '\n';
            ++input_pos;
            state_ = advance_state_on_newline(state_);
        }

        ring_memory_pool_ptr_->producer_set_finished();
    }

private:
    enum class State
    {
        ReadHeader,
        ReadSequence,
        ReadPlus,
        ReadQuality
    };

    State advance_state_on_newline(const State current) const
    {
        switch (current)
        {
        case State::ReadHeader:
            return State::ReadSequence;
        case State::ReadSequence:
            return State::ReadPlus;
        case State::ReadPlus:
            return State::ReadQuality;
        case State::ReadQuality:
            return State::ReadHeader;
        default:
            return State::ReadHeader;
        }
    }

    inline bool acquire_block(char *&block_ptr,
                              uint64_t &write_size,
                              uint64_t &last_newline_pos,
                              bool &has_block)
    {
        if (!ring_memory_pool_ptr_->producer_dequeue(block_ptr))
        {
            return false;
        }

        has_block = true;
        write_size = 0;
        last_newline_pos = kNoNewlineInBlock;
        assert(left_buffer_size_ <= sizeof(left_buffer_) && "left_buffer_size_ exceeds left_buffer_ capacity");

        if (left_buffer_size_ > 0)
        {
            std::memcpy(block_ptr, left_buffer_, left_buffer_size_);
            write_size = left_buffer_size_;
            left_buffer_size_ = 0;
        }
        return true;
    }

    inline void store_overlap_from_block_end(const char *block_ptr,
                                             const uint64_t write_size,
                                             const uint64_t last_newline_pos,
                                             const uint64_t overlap)
    {
        left_buffer_size_ = 0;
        if (overlap == 0 || write_size == 0)
        {
            return;
        }

        const uint64_t keep_window = (write_size < overlap) ? write_size : overlap;
        const uint64_t start = write_size - keep_window;

        if (last_newline_pos != kNoNewlineInBlock && last_newline_pos >= start)
        {
            const uint64_t keep = write_size - (last_newline_pos + 1);
            if (keep > 0)
            {
                std::memcpy(left_buffer_, block_ptr + last_newline_pos + 1, keep);
                left_buffer_size_ = static_cast<size_t>(keep);
            }
            return;
        }

        std::memcpy(left_buffer_, block_ptr + start, keep_window);
        left_buffer_size_ = static_cast<size_t>(keep_window);
    }

    inline void publish_current_block(char *&block_ptr,
                                      uint64_t &write_size,
                                      uint64_t &last_newline_pos,
                                      bool &has_block)
    {
        if (!has_block)
        {
            return;
        }
        if (write_size > 0)
        {
            ring_memory_pool_ptr_->producer_enqueue({block_ptr, write_size});
        }
        else
        {
            ring_memory_pool_ptr_->consumer_enqueue(block_ptr);
        }
        has_block = false;
        block_ptr = nullptr;
        write_size = 0;
        last_newline_pos = kNoNewlineInBlock;
    }

    static constexpr uint64_t kNoNewlineInBlock = static_cast<uint64_t>(-1);

    State state_ = State::ReadHeader;
    int fd_ = -1;
    off_t file_size_ = 0;
    off_t have_read_ = 0;
    uint64_t chunk_size_;
    std::string filename_;
    RingMemoryPool<RING_MEMORY_POOL_CAPACITY> *ring_memory_pool_ptr_;
    std::vector<char> file_buffer_; // 预留给读缓冲参数，保持接口兼容
    char left_buffer_[128];         // 块在碱基行中间截断时，保存最后(k-1)个碱基用于下块前缀
    size_t left_buffer_size_ = 0;   // left_buffer_中有效碱基字节数
};

#endif