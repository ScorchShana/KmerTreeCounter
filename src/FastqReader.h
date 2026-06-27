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
#include <cstdlib>
#include <zlib.h>

template <uint32_t N>
class FastqReader
{

    using content_type = std::pair<char*, uint64_t>;

    enum class State
    {
        ReadHeader,
        ReadSequence,
        ReadPlus,
        ReadQuality
    };

    static constexpr uint64_t kNoNewlineInBlock = static_cast<uint64_t>(-1);

    int acbs_index = 0;

    State state_ = State::ReadHeader;
    int fd_ = -1;
    off_t file_size_ = 0;
    off_t have_read_ = 0;
    uint64_t chunk_size_;
    std::string filename_;
    SPMCRingMemoryPool<READER_PARSER_RING_MEMORY_POOL_CAPACITY>* ring_memory_pool_ptr_;
    char* file_buffer;            // 预留给读缓冲参数，保持接口兼容
    char left_buffer_[128];       // 块在碱基行中间截断时，保存最后(k-1)个碱基用于下块前缀
    size_t left_buffer_size_ = 0; // left_buffer_中有效碱基字节数
    bool is_gz_file = false; // 是否为 gzip 压缩文件
    gzFile gzfile_ = nullptr;

public:
#ifdef TEST_MODE
    uint64_t total_dequeue_spin_time = 0;
    uint64_t total_enqueue_spin_time = 0;
    uint64_t aio_wait_spin_time = 0;
#endif

    int k;

    explicit FastqReader(const std::string& filename, const int in_k, uint64_t chunk_size, SPMCRingMemoryPool<READER_PARSER_RING_MEMORY_POOL_CAPACITY>* in_ring_memory_pool_ptr)
        : acbs_index(0), filename_(filename), ring_memory_pool_ptr_(in_ring_memory_pool_ptr),
        k(in_k), chunk_size_(chunk_size)
    {
        assert(ring_memory_pool_ptr_ != nullptr && "SPMCRingMemoryPool pointer must not be null");
        assert(k > 0 && "k must be positive");
        assert(k < 128 && "k must be < 128");
        assert(ring_memory_pool_ptr_->blockSize() > 1024 && "block_size must be > 1024");

        if (ring_memory_pool_ptr_ == nullptr)
        {
            std::cerr << "RingMemoryPool pointer is null" << std::endl;
            std::exit(-1);
        }

        fd_ = ::open(filename.data(), O_RDONLY);
        if (fd_ == -1)
        {
            std::cerr << "Failed to open file" << std::endl;
            std::exit(-1);
        }
        struct stat st;
        if (fstat(fd_, &st) == -1)
        {
            std::cerr << "Failed to get file size" << std::endl;
            std::exit(-1);
        }
        file_size_ = st.st_size;

        //判断是否为gz文件
        unsigned char buf[2];
        ssize_t n = ::read(fd_, buf, 2);
        is_gz_file = (n == 2 && buf[0] == 0x1F && buf[1] == 0x8B);

        if(is_gz_file)
        {
            gzfile_ = gzopen(filename.c_str(), "rb");
            if (gzfile_ == nullptr) {
                std::cerr << "Failed to open gzip file: " << filename << std::endl;
                std::exit(-1);
            }
        }

        else
        {
            if (lseek(fd_, 0, SEEK_SET) == -1) {
                std::cerr << "Failed to seek to beginning" << std::endl;
                std::exit(-1);
            }

            posix_fadvise(fd_, 0, 0, POSIX_FADV_SEQUENTIAL); // 顺序访问提示
        }

        file_buffer = static_cast<char*>(std::aligned_alloc(4096, chunk_size_));
        if (!file_buffer)
        {
            std::cerr << "Failed to allocate aligned memory for file buffe9r: " << std::endl;
            std::exit(-1);
        }
    }

    ~FastqReader()
    {
        if (gzfile_ != nullptr) 
        {
            gzclose(gzfile_);
        }
            if (fd_ != -1)
        {
            ::close(fd_);
        }
        if (file_buffer)
        {
            std::free(file_buffer);
        }
    }

    void read()
    {
        assert(ring_memory_pool_ptr_ != nullptr && "RingMemoryPool pointer must not be null");
        char* block_ptr = nullptr;
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

        double cur_percent = 0.0;
        int last_reported_percent = 0;


        while (!stop)
        {
            if (!has_block)
            {
                acquire_block(block_ptr, write_size, last_newline_pos, has_block);
            }

            if (input_pos >= input_size && !eof)
            {
                ssize_t bytes_read = 0;
                if(is_gz_file){
                    bytes_read = gzread(gzfile_, file_buffer, chunk_size_);
                }else{
                    bytes_read = ::read(fd_, file_buffer, chunk_size_);
                }
                
                if (bytes_read < 0) [[unlikely]]
                {
                    std::cerr << "Failed to read fastq data" << std::endl;
                    std::exit(-1);
                }

                if (bytes_read == 0) [[unlikely]]
                {
                    eof = true;
                }
                else
                {
                    input_pos = 0;
                    input_size = static_cast<uint64_t>(bytes_read);
                    have_read_ += bytes_read;
               
                    // 进度更新频率为每增加1%，或达到100%时更新一次
                    double current_read = is_gz_file ? static_cast<double>(gzoffset(gzfile_)) : static_cast<double>(have_read_);
                    cur_percent = current_read / static_cast<double>(file_size_);
                    int cur_percent_int = static_cast<int>(cur_percent * 100);
                    if (cur_percent_int - last_reported_percent >= 1 || cur_percent_int == 100)
                    {
                        std::cout << "\rFastqReader progress: " << cur_percent_int << "%";
                    }
                    if (cur_percent_int >= 100) [[unlikely]]
                    {
                        std::cout << std::endl;
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

            const char* input_begin = file_buffer;

            if (state_ != State::ReadSequence)
            {
                const char* cur = input_begin + input_pos;
                const uint64_t remain = input_size - input_pos;
                const void* newline_ptr = std::memchr(cur, '\n', remain);

                if (newline_ptr == nullptr)
                {
                    input_pos = input_size;
                    continue;
                }

                const char* newline = static_cast<const char*>(newline_ptr);
                input_pos = static_cast<uint64_t>(newline - input_begin) + 1;
                state_ = advance_state_on_newline(state_);
                continue;
            }

            // HOT PATH - sequence line is copied in batches instead of char-by-char.
            const char* seq_cur = input_begin + input_pos;
            const uint64_t seq_remain = input_size - input_pos;
            const void* newline_ptr = std::memchr(seq_cur, '\n', seq_remain);
            const uint64_t sequence_len = (newline_ptr == nullptr)
                ? seq_remain
                : static_cast<uint64_t>(static_cast<const char*>(newline_ptr) - seq_cur);

            uint64_t copied = 0;
            while (copied < sequence_len)
            {
                if (write_size == block_size)
                {
                    store_overlap_from_block_end(block_ptr, write_size, last_newline_pos, overlap);
                    publish_current_block(block_ptr, write_size, last_newline_pos, has_block);
                    acquire_block(block_ptr, write_size, last_newline_pos, has_block);
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
                acquire_block(block_ptr, write_size, last_newline_pos, has_block);
            }

            last_newline_pos = write_size;
            block_ptr[write_size++] = '\n';
            ++input_pos;
            state_ = advance_state_on_newline(state_);
        }

        ring_memory_pool_ptr_->producer_set_finished();
    }

private:
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

    void acquire_block(char*& block_ptr,
        uint64_t& write_size,
        uint64_t& last_newline_pos,
        bool& has_block)
    {

        // int spin_count = 1;
        // int backoff_iterations = 1;

        // while (!ring_memory_pool_ptr_->producer_try_dequeue(block_ptr))
        // {
        //     spin_time++;
        // }
        ring_memory_pool_ptr_->producer_dequeue(block_ptr);

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
    }

    inline void store_overlap_from_block_end(const char* block_ptr,
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

    inline void publish_current_block(char*& block_ptr,
        uint64_t& write_size,
        uint64_t& last_newline_pos,
        bool& has_block)
    {
        if (!has_block)
        {
            return;
        }
        if (write_size > 0)
        {
            ring_memory_pool_ptr_->producer_enqueue({ block_ptr, write_size });
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
};

#endif