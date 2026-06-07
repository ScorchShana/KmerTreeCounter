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

/*template <uint32_t N>
class alignas(PAGE_SIZE) ExportWriter
{
    constexpr static int MAX_SPIN_TIME = 1024;
    constexpr static int MAX_BACKOFF = 128;

    static constexpr uint64_t AIO_QUEUE_DEPTH = 128;
    static constexpr uint64_t PAGE_POOL_SIZE = EXPORT_FILES_SIZE + 2 * AIO_QUEUE_DEPTH;
    static constexpr uint64_t PER_PAGE_STORAGE_SIZE = 256 * 1024ULL;

private:
    // 连接到 FastqParser，用于获取含有低频 k-mer 的数据块
    RingMemoryPool<EXPORT_RING_MEMORY_POOL_CAPACITY> *ring_memory_pool = nullptr;

    //页存储区，直接内嵌，零额外堆分配 ===
    alignas(PAGE_SIZE) char page_storage[PAGE_POOL_SIZE * PER_PAGE_STORAGE_SIZE];

    // 空闲页栈：存储可用页在 page_storage 中的索引
    std::array<uint16_t, PAGE_POOL_SIZE> free_stack{};
    std::atomic<uint64_t> free_top{PAGE_POOL_SIZE};

    // 每个文件的当前写入页（可能为 nullptr，表示尚未分配）
    std::array<char *, EXPORT_FILES_SIZE> cur_pages{};
    std::array<uint64_t, EXPORT_FILES_SIZE> cur_used{};
    std::array<uint64_t, EXPORT_FILES_SIZE> file_offsets{};
    std::array<uint64_t, EXPORT_FILES_SIZE> total_counts{};
    std::array<int, EXPORT_FILES_SIZE> fds{};

    io_context_t aio_ctx = 0;

    // 专门执行落盘操作的单线程
    std::thread worker_thread;
    std::thread reaper_thread;
    std::atomic<bool> is_running{false};
    std::atomic<bool> reaper_running{false};

public:
    ExportWriter() = delete;
    ExportWriter(const ExportWriter &) = delete;
    ExportWriter &operator=(const ExportWriter &) = delete;
    ExportWriter(ExportWriter &&) = delete;
    ExportWriter &operator=(ExportWriter &&) = delete;

    explicit ExportWriter(RingMemoryPool<EXPORT_RING_MEMORY_POOL_CAPACITY> *pool)
        : ring_memory_pool(pool)
    {
        // 初始化页池：所有页索引入栈，并清零
        for (uint64_t i = 0; i < PAGE_POOL_SIZE; ++i)
        {
            free_stack[i] = static_cast<uint16_t>(i);
            // std::memset(page_storage + i * PER_PAGE_STORAGE_SIZE, 0, PER_PAGE_STORAGE_SIZE);
        }
        free_top = PAGE_POOL_SIZE;

        // 初始化 libaio
        if (io_setup(AIO_QUEUE_DEPTH, &aio_ctx) != 0)
        {
            std::cerr << "Failed to setup aio context" << std::endl;
            std::exit(-1);
        }

        // 打开文件，写入 4KB 空 header
        for (uint64_t i = 0; i < EXPORT_FILES_SIZE; ++i)
        {
            const std::string filename = temp_dir + "low_" + std::to_string(i) + ".bin";
            fds[i] = ::open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_DIRECT, 0644);
            if (fds[i] < 0) [[unlikely]]
            {
                std::cerr << "Failed to open export file: " << filename << std::endl;
                std::exit(-1);
            }

            char *header_page = alloc_page();
            ssize_t ret = ::pwrite(fds[i], header_page, PER_PAGE_STORAGE_SIZE, 0);
            if (ret != static_cast<ssize_t>(PER_PAGE_STORAGE_SIZE)) [[unlikely]]
            {
                std::cerr << "Failed to write header for file: " << filename << std::endl;
                std::exit(-1);
            }
            release_page(header_page);

            file_offsets[i] = PER_PAGE_STORAGE_SIZE;
            cur_pages[i] = nullptr;
            cur_used[i] = 0;
            total_counts[i] = 0;

        }
    }

    ~ExportWriter()
    {
        for (uint64_t i = 0; i < EXPORT_FILES_SIZE; ++i)
        {
            if (fds[i] >= 0)
            {
                ::close(fds[i]);
                fds[i] = -1;
            }
        }
        if (aio_ctx != 0)
        {
            io_destroy(aio_ctx);
            aio_ctx = 0;
        }
    }

    void start()
    {
        reaper_running = true;
        is_running = true;
        reaper_thread = std::thread(&ExportWriter::reaper_loop, this);
        worker_thread = std::thread(&ExportWriter::writer_loop, this);
    }

    void join()
    {
        // 等待写线程结束（writer_loop 内部会调用 flush_all_buffers）
        if (is_running && worker_thread.joinable())
        {
            worker_thread.join();
            is_running = false;
        }

        // 刷空所有剩余 page 并等待 AIO 完成（收割线程仍运行中）
        flush_all_buffers();

        // 停收割线程（退出前会 drain 完所有遗留 AIO）
        reaper_running = false;
        if (reaper_thread.joinable())
            reaper_thread.join();

        write_headers();

        for (uint64_t i = 0; i < EXPORT_FILES_SIZE; ++i)
        {
            if (fds[i] >= 0)
            {
                ::fsync(fds[i]);
                ::close(fds[i]);
                fds[i] = -1;
            }
        }

        if (aio_ctx != 0)
        {
            io_destroy(aio_ctx);
            aio_ctx = 0;
        }
    }

private:
    char *alloc_page()
    {
        while (free_top.load(std::memory_order_acquire) == 0)
        {
            cpu_relax();
        }
        uint64_t idx = free_stack[free_top.fetch_sub(1, std::memory_order_acq_rel) - 1];
        return page_storage + idx * PER_PAGE_STORAGE_SIZE;
    }

    void release_page(char *p)
    {
        uint64_t idx = static_cast<uint64_t>(p - page_storage) / PER_PAGE_STORAGE_SIZE;
        uint64_t top = free_top.fetch_add(1, std::memory_order_release);
        free_stack[top] = static_cast<uint16_t>(idx);
    }

    void reaper_loop()
    {
        struct io_event events[AIO_QUEUE_DEPTH];
        struct timespec timeout = {0, 50 * 1000000};

        while (reaper_running.load(std::memory_order_acquire))
        {
            int n = io_getevents(aio_ctx, 1, static_cast<long>(AIO_QUEUE_DEPTH), events, &timeout);
            for (int i = 0; i < n; ++i)
            {
                if (events[i].res < 0) [[unlikely]]
                {
                    std::cerr << "AIO write failed: " << events[i].res << std::endl;
                    std::exit(-1);
                }
                release_page(static_cast<char *>(events[i].data));
            }
        }

        struct io_event drain_events[AIO_QUEUE_DEPTH];
        while (true)
        {
            int n = io_getevents(aio_ctx, 0, static_cast<long>(AIO_QUEUE_DEPTH), drain_events, nullptr);
            if (n == 0) break;
            for (int i = 0; i < n; ++i)
                release_page(static_cast<char *>(drain_events[i].data));
        }
    }

    void submit_page(uint64_t prefix)
    {
        char *page = cur_pages[prefix];
        if (page == nullptr)
            return;

        struct iocb cb;
        io_prep_pwrite(&cb, fds[prefix], page, PER_PAGE_STORAGE_SIZE, static_cast<long long>(file_offsets[prefix]));
        cb.data = page;

        struct iocb *cbs[1] = {&cb};
        while (true)
        {
            int ret = io_submit(aio_ctx, 1, cbs);
            if (ret == 1)
                break;
            if (ret < 0 && ret == -EAGAIN)
            {
                cpu_relax();
                continue;
            }
            std::cerr << "io_submit failed: " << ret << std::endl;
            std::exit(-1);
        }

        file_offsets[prefix] += PER_PAGE_STORAGE_SIZE;
        cur_pages[prefix] = nullptr;
        cur_used[prefix] = 0;
    }

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
            const uint64_t file_index = get_file_index(val);
            total_counts[file_index]++;

            char *&page = cur_pages[file_index];
            uint64_t &used = cur_used[file_index];

            if (page == nullptr)
            {
                page = alloc_page();
                used = 0;
            }

            if (used + sizeof(kmer<N>) > PER_PAGE_STORAGE_SIZE)
            {
                submit_page(file_index);
                page = alloc_page();
                used = 0;
            }

            std::memcpy(page + used, &val, sizeof(kmer<N>));
            used += sizeof(kmer<N>);
        }
    }

    uint64_t get_file_index(const kmer<N> &val) const
    {
        constexpr uint32_t shift = 2 * (BASES_PER_U64T - ROOT_BASES);
        return val.data[0] >> shift;
    }

    void flush_all_buffers()
    {
        for (uint64_t i = 0; i < EXPORT_FILES_SIZE; ++i)
        {
            if (cur_used[i] > 0)
            {
                submit_page(i);
            }
            else if (cur_pages[i] != nullptr)
            {
                release_page(cur_pages[i]);
                cur_pages[i] = nullptr;
                cur_used[i] = 0;
            }
        }
        // 收割线程会持续调用 io_getevents 归还页面
        while (free_top.load(std::memory_order_acquire) < PAGE_POOL_SIZE)
        {
            cpu_relax();
        }
    }

    void write_headers()
    {
        for (uint64_t i = 0; i < EXPORT_FILES_SIZE; ++i)
        {
            char *page = alloc_page();
            std::memset(page, 0, PAGE_SIZE);
            *reinterpret_cast<uint64_t *>(page) = total_counts[i];
            ssize_t ret = ::pwrite(fds[i], page, PAGE_SIZE, 0);
            if (ret != static_cast<ssize_t>(PAGE_SIZE)) [[unlikely]]
            {
                std::cerr << "Failed to write header back for file " << i << std::endl;
                std::exit(-1);
            }
            release_page(page);
        }
    }
};
*/

/*
template <uint32_t N>
class ExportWriter
{

    constexpr static int MAX_SPIN_TIME = 1024;
    constexpr static int MAX_BACKOFF = 128;

private:
    // 连接到 FastqParser，用于获取含有低频 k-mer 的数据块
    RingMemoryPool<EXPORT_RING_MEMORY_POOL_CAPACITY> *ring_memory_pool = nullptr;
    // 为每个根节点（或者短前缀）维护一个文件句柄
    std::array<int, EXPORT_FILES_SIZE> files{};
    // 每个根节点对应的最大缓冲容量（默认 512KB）
    static constexpr uint64_t ROOT_BUFFER_KMER_CAPACITY = EXPORT_ROOT_BUFFER_SIZE / sizeof(kmer<N>);
    static constexpr uint64_t LOCAL_BUFFERS_CAPACITY = EXPORT_RING_MEMORY_POOL_BLOCK_SIZE / sizeof(kmer<N>);
    // 对应每个根节点的本地缓存，当装满时通过 fwrite 顺序落盘
    std::array<kmer<N> *, EXPORT_FILES_SIZE> root_buffers{};
    std::array<uint32_t, EXPORT_FILES_SIZE> root_buffer_counts{};

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
            const std::string filename = temp_dir + "low_" + std::to_string(i) + ".bin";
            files[i] = ::open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
            if (files[i] < 0) [[unlikely]]
            {
                std::cerr << "Failed to open export file: " << filename << std::endl;
                std::exit(-1);
            }
        }
        std::memset(root_buffer_counts.data(), 0, root_buffer_counts.size() * sizeof(uint32_t));
    }

    ~ExportWriter()
    {
        for (uint64_t i = 0; i < files.size(); ++i)
        {
            if (files[i] >= 0)
            {
                ::close(files[i]);
                files[i] = -1;
            }
            delete[] root_buffers[i];
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
            if (files[i] >= 0)
            {
                ::close(files[i]);
                files[i] = -1;
            }
        }
    }

private:
    void writer_loop()
    {

        for (uint64_t i = 0; i < files.size(); ++i)
        {
            root_buffers[i] = new kmer<N>[ROOT_BUFFER_KMER_CAPACITY];
        }

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
            auto &buffer_count = root_buffer_counts[root_prefix];
            buffer[buffer_count++] = val;
            if (buffer_count >= ROOT_BUFFER_KMER_CAPACITY) [[unlikely]]
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
        auto &buffer_count = root_buffer_counts[root_prefix];

        if (::write(files[root_prefix], buffer, sizeof(kmer<N>) * buffer_count) != sizeof(kmer<N>) * buffer_count) [[unlikely]]
        {
            std::cerr << "Failed to write k-mer data" << std::endl;
            std::exit(-1);
        }
        buffer_count = 0;
    }

    void flush_all_buffers()
    {
        for (uint64_t i = 0; i < root_buffers.size(); ++i)
        {
            if (root_buffer_counts[i] != 0)
            {
                flush_root_buffer(i);
            }
        }
    }
};
*/

template <uint32_t N>
class ExportWriter
{

    constexpr static int MAX_SPIN_TIME = 1024;
    constexpr static int MAX_BACKOFF = 128;

    constexpr static uint64_t RAW_BUFFER_SIZE = 2ULL * 1024 * 1024; // 2MB 原始数据块大小
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
        if (::write(fd, buffer, sizeof(kmer<N>) * buffer_count) != sizeof(kmer<N>) * buffer_count) [[unlikely]]
        {
            std::cerr << "Failed to write k-mer data" << std::endl;
            std::exit(-1);
        }
        buffer_count = 0;
    }
};

#endif // EXPORT_WRITER_HEADER
