#ifndef CONCURRENT_MEMORY_POOL_HEADER
#define CONCURRENT_MEMORY_POOL_HEADER

// 并发内存池 - NUMA感知设计
// 固定块模型（4KB），线程本地缓存，跨线程释放支持

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <mutex>
#include <thread>
#include <vector>
#include <stdexcept>
#include <new>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include "definition.h"

// 系统头文件
#include <sys/mman.h>
#include <sched.h>
#include <unistd.h>
#include <errno.h>

// NUMA 支持 - 自动检测 libnuma 是否可用
#if __has_include(<numa.h>)
#include <numa.h>
#include <numaif.h>
#define HAS_LIBNUMA 1
#else
#define HAS_LIBNUMA 0
#endif

//==============================================================================
// 可配置常量
//==============================================================================

// 最大 NUMA 节点数
static constexpr int MAX_NUMA_NODES = 16;

// 固定块大小（4KB）
static constexpr size_t BLOCK_SIZE = 4096;

// 从 Arena 批量获取的块数
static constexpr size_t BATCH_SIZE = 64;

// 线程本地栈最大缓存块数
static constexpr size_t TLS_CAPACITY = 128;

// 每个远程链表的最大长度
static constexpr size_t REMOTE_LIST_CAPACITY = 32;

//static constexpr size_t CACHE_LINE_SIZE = 64;



//==============================================================================
// 辅助工具
//==============================================================================

// 向上对齐到指定边界
constexpr size_t align_up(size_t value, size_t alignment)
{
    return (value + alignment - 1) & ~(alignment - 1);
}

// 检查是否对齐
constexpr bool is_aligned(size_t value, size_t alignment)
{
    return (value & (alignment - 1)) == 0;
}

//==============================================================================
// FreeBlock - 侵入式链表节点
//==============================================================================

// 空闲块头部结构，用于链接空闲块
// 当块被分配时，此空间可供用户使用
struct FreeBlock
{
    FreeBlock *next;
};

//==============================================================================
// Arena - NUMA 节点内存管理
//==============================================================================

// 每个 NUMA 节点一个 Arena，管理本节点的内存块
// 使用 alignas 避免伪共享
struct alignas(CACHE_LINE_SIZE) Arena
{
    void *start_addr;             // 本 Arena 的起始地址
    void *end_addr;               // 本 Arena 的结束地址（不包含）
    std::atomic<char *> bump_cursor;
    FreeBlock *central_free_list; // 中央空闲链表
    std::mutex mutex;             // 保护中央空闲链表的互斥锁

    Arena() : start_addr(nullptr), end_addr(nullptr), bump_cursor(nullptr), central_free_list(nullptr) {}

    // 检查指针是否属于此 Arena
    bool contains(void *ptr) const
    {
        return ptr >= start_addr && ptr < end_addr;
    }

    // 获取此 Arena 管理的块数
    size_t block_count() const
    {
        return static_cast<size_t>(
                   static_cast<char *>(end_addr) - static_cast<char *>(start_addr)) /
               BLOCK_SIZE;
    }
};

//==============================================================================
// RemoteListHead - 远程链表头
//==============================================================================

// 用于 ThreadLocalCache 中的远程链表
struct RemoteListHead
{
    FreeBlock *head;
    size_t count;

    RemoteListHead() : head(nullptr), count(0) {}
};

//==============================================================================
// ThreadLocalCache - 线程本地缓存
//==============================================================================

struct ThreadLocalCache
{
    // 本地空闲栈（使用侵入式链表实现栈）
    FreeBlock *local_free_stack;
    size_t local_free_count;

    // 远程链表数组，每个 Arena 一个
    // 用于暂存属于其他 Arena 的块
    alignas(CACHE_LINE_SIZE) RemoteListHead remote_lists[MAX_NUMA_NODES];

    // 指向本地 Arena 的指针
    Arena *local_arena;
    int local_arena_index;

    // 所属的 ConcurrentMemoryPool 指针（用于归还内存）
    class ConcurrentMemoryPool *pool;

    ThreadLocalCache()
        : local_free_stack(nullptr), local_free_count(0), local_arena(nullptr), local_arena_index(-1), pool(nullptr) {}

    // 析构时归还所有缓存
    ~ThreadLocalCache();
};

//==============================================================================
// MemoryPool - 顶层内存池类
//==============================================================================

class ConcurrentMemoryPool
{
public:
    // 构造函数
    // total_bytes: 总字节数（内部会按 BLOCK_SIZE 对齐）
    // num_threads: 初始化 first-touch 的线程数
    explicit ConcurrentMemoryPool(size_t total_bytes, size_t num_threads);

    // 禁止拷贝与赋值
    ConcurrentMemoryPool(const ConcurrentMemoryPool &) = delete;
    ConcurrentMemoryPool &operator=(const ConcurrentMemoryPool &) = delete;

    // 析构函数，释放所有 mmap 内存
    ~ConcurrentMemoryPool();

    // 分配一个 4KB 块，失败时抛出 std::bad_alloc
    void *allocate();

    // 分配可变大小的块，失败时抛出 std::bad_alloc
    void *allocate_large(size_t bytes);

    // 释放一个之前分配的块
    void deallocate(void *ptr);

    // 获取当前线程所属 Arena 的索引
    int get_current_arena_index() const;

    // 将当前线程的 TLS 缓存中的所有块归还给各自的 Arena
    static void reclaim_thread_cache();

    // 获取线程本地缓存
    static ThreadLocalCache &get_thread_local_cache();

private:
    //--------------------------------------------------------------------------
    // 内部辅助函数
    //--------------------------------------------------------------------------

    // 初始化 NUMA 信息
    void init_numa_info();

    // 分配内存（尝试大页）并设置 Arena 地址范围
    void allocate_memory(size_t total_bytes);

    // 执行 first-touch 预加载
    void perform_first_touch(size_t num_threads);

    // 从 Arena 批量获取块，返回实际获取块数（0 表示获取失败）
    size_t batch_allocate_from_arena(Arena &arena, FreeBlock **out_head, FreeBlock **out_tail, size_t count);

    char *bump_allocate_from_arena(Arena &arena, size_t bytes);

    FreeBlock *batch_allocate_from_bump(Arena &arena, size_t &out_count, FreeBlock **out_tail);

    // 批量归还块到 Arena
    void batch_deallocate_to_arena(Arena &arena, FreeBlock *head, FreeBlock *tail, size_t count);

    // 根据地址查找 Arena 索引
    int find_arena_index(void *ptr) const;

    // 获取当前线程应该使用的 Arena 索引
    int get_thread_arena_index() const;

    //--------------------------------------------------------------------------
    // 成员变量
    //--------------------------------------------------------------------------

    // Arena 数组
    alignas(CACHE_LINE_SIZE) Arena arenas_[MAX_NUMA_NODES];

    // 实际使用的 Arena 数量
    int num_arenas_;

    // 内存区域起始地址
    void *memory_start_;

    // 内存区域总大小
    size_t memory_size_;

    // 实际使用的页大小（4KB 或 2MB）
    size_t page_size_;

    // 总块数
    size_t total_blocks_;

    // NUMA 是否可用
    bool numa_available_;

    // 每个 NUMA 节点的 CPU 数量
    int cpus_per_node_[MAX_NUMA_NODES];

    // 总可用 CPU 数
    int total_cpus_;

    // 线程本地缓存
    static inline thread_local ThreadLocalCache tls_cache_;

    // 友元类
    friend struct ThreadLocalCache;
};

//==============================================================================
// ThreadLocalCache 析构函数实现
//==============================================================================

inline ThreadLocalCache::~ThreadLocalCache()
{
    if (!pool)
        return;

    // 归还本地空闲栈中的所有块
    while (local_free_stack)
    {
        FreeBlock *block = local_free_stack;
        local_free_stack = block->next;
        local_free_count--;

        // 直接归还到对应的 Arena
        int arena_idx = pool->find_arena_index(block);
        if (arena_idx >= 0 && arena_idx < pool->num_arenas_)
        {
            Arena &arena = pool->arenas_[arena_idx];
            std::lock_guard<std::mutex> lock(arena.mutex);
            block->next = arena.central_free_list;
            arena.central_free_list = block;
        }
    }

    // 归还所有远程链表中的块
    for (int i = 0; i < MAX_NUMA_NODES; ++i)
    {
        RemoteListHead &rl = remote_lists[i];
        if (rl.head)
        {
            if (i < pool->num_arenas_)
            {
                Arena &arena = pool->arenas_[i];
                std::lock_guard<std::mutex> lock(arena.mutex);
                // 将整个链表挂到中央空闲链表
                FreeBlock *tail = rl.head;
                while (tail->next)
                {
                    tail = tail->next;
                }
                tail->next = arena.central_free_list;
                arena.central_free_list = rl.head;
            }
            rl.head = nullptr;
            rl.count = 0;
        }
    }
}

//==============================================================================
// MemoryPool 实现
//==============================================================================

inline ConcurrentMemoryPool::ConcurrentMemoryPool(size_t total_bytes, size_t num_threads)
    : num_arenas_(1), memory_start_(nullptr), memory_size_(0), page_size_(4096), total_blocks_(0), numa_available_(false), total_cpus_(1)
{
    // 初始化数组
    std::memset(cpus_per_node_, 0, sizeof(cpus_per_node_));

    // 初始化 NUMA 信息
    init_numa_info();

    // 计算总块数（按 BLOCK_SIZE 对齐）
    total_blocks_ = total_bytes / BLOCK_SIZE;
    if (total_blocks_ == 0)
    {
        std::cerr << "total_bytes must be at least BLOCK_SIZE (4096)" << std::endl;
        std::exit(-1);
    }

    // 分配内存并设置 Arena 地址范围
    allocate_memory(total_blocks_ * BLOCK_SIZE);

    // 执行 first-touch 预加载
    perform_first_touch(num_threads);
}

inline ConcurrentMemoryPool::~ConcurrentMemoryPool()
{
    // 清理当前线程的 TLS 缓存（如果有的话）
    // 注意：TLS 析构会在 ConcurrentMemoryPool 析构之后发生，所以需要先清理
    ThreadLocalCache &tls = tls_cache_;
    if (tls.pool == this)
    {
        // 清空本地空闲栈，不归还到 Arena（因为 Arena 即将销毁）
        tls.local_free_stack = nullptr;
        tls.local_free_count = 0;
        tls.local_arena = nullptr;
        tls.local_arena_index = -1;

        // 清空远程链表
        for (int i = 0; i < MAX_NUMA_NODES; ++i)
        {
            tls.remote_lists[i].head = nullptr;
            tls.remote_lists[i].count = 0;
        }

        tls.pool = nullptr;
    }

    if (memory_start_)
    {
        munmap(memory_start_, memory_size_);
        memory_start_ = nullptr;
    }
}

inline void ConcurrentMemoryPool::init_numa_info()
{
#if HAS_LIBNUMA
    if (numa_available() == -1)
    {
        numa_available_ = false;
        num_arenas_ = 1;
        cpus_per_node_[0] = 1;
        total_cpus_ = 1;
        return;
    }

    numa_available_ = true;
    num_arenas_ = numa_num_configured_nodes();
    if (num_arenas_ > MAX_NUMA_NODES)
    {
        num_arenas_ = MAX_NUMA_NODES;
    }

    // 获取当前进程的 CPU 亲和性
    cpu_set_t mask;
    CPU_ZERO(&mask);
    if (sched_getaffinity(0, sizeof(mask), &mask) == 0)
    {
        // 统计每个 NUMA 节点的可用 CPU 数
        total_cpus_ = 0;
        for (int node = 0; node < num_arenas_; ++node)
        {
            int count = 0;
            for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
            {
                if (CPU_ISSET(cpu, &mask))
                {
                    if (numa_node_of_cpu(cpu) == node)
                    {
                        ++count;
                    }
                }
            }
            cpus_per_node_[node] = count;
            total_cpus_ += count;
        }

        // 如果没有找到 CPU，使用默认值
        if (total_cpus_ == 0)
        {
            total_cpus_ = 1;
            cpus_per_node_[0] = 1;
        }
    }
    else
    {
        // 无法获取亲和性，使用默认值
        total_cpus_ = 1;
        cpus_per_node_[0] = 1;
    }
#else
    numa_available_ = false;
    num_arenas_ = 1;
    cpus_per_node_[0] = 1;
    total_cpus_ = 1;
#endif
}

inline void ConcurrentMemoryPool::allocate_memory(size_t total_bytes)
{
    // 确保总大小是页大小的整数倍
    // 首先尝试 2MB 大页
    size_t huge_page_size = 2 * 1024 * 1024; // 2MB

    // 计算需要的总大小（对齐到页大小）
    size_t aligned_size = align_up(total_bytes, huge_page_size);

    // 尝试使用大页
    void *addr = mmap(nullptr, aligned_size,
                      PROT_READ | PROT_WRITE,
                      MAP_PRIVATE | MAP_ANONYMOUS | MAP_HUGETLB,
                      -1, 0);

    if (addr != MAP_FAILED)
    {
        // 大页分配成功
        memory_start_ = addr;
        memory_size_ = aligned_size;
        page_size_ = huge_page_size;
    }
    else
    {
        // 大页分配失败，回退到普通页
        aligned_size = align_up(total_bytes, static_cast<size_t>(sysconf(_SC_PAGESIZE)));
        addr = mmap(nullptr, aligned_size,
                    PROT_READ | PROT_WRITE,
                    MAP_PRIVATE | MAP_ANONYMOUS,
                    -1, 0);

        if (addr == MAP_FAILED)
        {
            std::cerr << "std::bad_alloc" << std::endl;
            std::exit(-1);
        }

        memory_start_ = addr;
        memory_size_ = aligned_size;
        page_size_ = static_cast<size_t>(sysconf(_SC_PAGESIZE));

        // 尝试建议使用透明大页
        madvise(memory_start_, memory_size_, MADV_HUGEPAGE);
    }

    // 设置每个 Arena 的地址范围
    char *current_addr = static_cast<char *>(memory_start_);
    size_t remaining_blocks = total_blocks_;

    for (int i = 0; i < num_arenas_; ++i)
    {
        size_t arena_blocks;
        if (i == num_arenas_ - 1)
        {
            arena_blocks = remaining_blocks;
        }
        else
        {
            arena_blocks = (total_blocks_ * cpus_per_node_[i]) / total_cpus_;
            if (arena_blocks > remaining_blocks)
            {
                arena_blocks = remaining_blocks;
            }
        }

        // 确保 Arena 大小是页大小的整数倍
        size_t arena_bytes = arena_blocks * BLOCK_SIZE;
        size_t aligned_arena_bytes = align_up(arena_bytes, page_size_);

        arenas_[i].start_addr = current_addr;
        arenas_[i].end_addr = current_addr + arena_bytes;
        arenas_[i].bump_cursor.store(current_addr, std::memory_order_relaxed);

        current_addr += aligned_arena_bytes;
        remaining_blocks -= arena_blocks;
    }
}

inline void ConcurrentMemoryPool::perform_first_touch(size_t num_threads)
{
    if (num_threads == 0)
    {
        num_threads = 1;
    }

    // 计算每个 Arena 分配的线程数
    std::vector<int> threads_per_arena(num_arenas_);
    int remaining_threads = static_cast<int>(num_threads);

    for (int i = 0; i < num_arenas_; ++i)
    {
        if (i == num_arenas_ - 1)
        {
            threads_per_arena[i] = remaining_threads;
        }
        else
        {
            threads_per_arena[i] = (num_threads * cpus_per_node_[i]) / total_cpus_;
            if (threads_per_arena[i] < 1)
                threads_per_arena[i] = 1;
            remaining_threads -= threads_per_arena[i];
        }
    }

    // 创建线程执行 first-touch
    std::vector<std::thread> threads;

    for (int arena_idx = 0; arena_idx < num_arenas_; ++arena_idx)
    {
        for (int rank = 0; rank < threads_per_arena[arena_idx]; ++rank)
        {
            threads.emplace_back([this, arena_idx, rank, &threads_per_arena]()
                                 {
            //--------------------------------------------------------------------------
            // 设置 CPU 亲和性：将线程绑定到对应的 NUMA 节点
            //--------------------------------------------------------------------------
#if HAS_LIBNUMA
                if (numa_available_) {
                    // 使用 numa_run_on_node 将线程绑定到指定 NUMA 节点
                    // 这确保 first-touch 发生在正确的节点上
                    numa_run_on_node(arena_idx);
                }
#endif
                
                Arena& arena = arenas_[arena_idx];
                
                // 计算此 Arena 的页数
                size_t arena_bytes = static_cast<char*>(arena.end_addr) - 
                                     static_cast<char*>(arena.start_addr);
                size_t pages_in_arena = arena_bytes / page_size_;
                if (pages_in_arena == 0) pages_in_arena = 1;
                
                // 计算此线程负责的页范围
                size_t pages_per_thread = pages_in_arena / threads_per_arena[arena_idx];
                if (pages_per_thread == 0) pages_per_thread = 1;
                
                size_t start_page = rank * pages_per_thread;
                size_t end_page = (rank == threads_per_arena[arena_idx] - 1) 
                                  ? pages_in_arena 
                                  : (rank + 1) * pages_per_thread;
                
                // 对每一页执行 first-touch
                char* page_start = static_cast<char*>(arena.start_addr) + start_page * page_size_;
                char* arena_end = static_cast<char*>(arena.end_addr);
                
                for (size_t page = start_page; page < end_page && page_start < arena_end; ++page) {
                    // First-touch: 写入一个字节
                    *page_start = 0;
                    page_start += page_size_;
                }
                });
        }
    }

    // 等待所有线程完成
    for (auto &t : threads)
    {
        t.join();
    }

}

inline char *ConcurrentMemoryPool::bump_allocate_from_arena(Arena &arena, size_t bytes)
{
    if (bytes == 0)
    {
        return nullptr;
    }

    size_t aligned_bytes = align_up(bytes, BLOCK_SIZE);
    char *arena_end = static_cast<char *>(arena.end_addr);

    for (;;)
    {
        char *cursor = arena.bump_cursor.load(std::memory_order_relaxed);
        if (!cursor)
        {
            return nullptr;
        }

        char *next = cursor + aligned_bytes;
        if (next > arena_end)
        {
            return nullptr;
        }

        if (arena.bump_cursor.compare_exchange_weak(cursor, next,
                                                    std::memory_order_relaxed,
                                                    std::memory_order_relaxed))
        {
            return cursor;
        }
    }
}

inline FreeBlock *ConcurrentMemoryPool::batch_allocate_from_bump(Arena &arena, size_t &out_count, FreeBlock **out_tail)
{
    if (out_count == 0)
    {
        return nullptr;
    }

    char *arena_end = static_cast<char *>(arena.end_addr);

    for (;;)
    {
        char *cursor = arena.bump_cursor.load(std::memory_order_relaxed);
        if (!cursor || cursor >= arena_end)
        {
            out_count = 0;
            return nullptr;
        }

        size_t available_bytes = static_cast<size_t>(arena_end - cursor);
        size_t available_blocks = available_bytes / BLOCK_SIZE;
        if (available_blocks == 0)
        {
            out_count = 0;
            return nullptr;
        }

        size_t allocate_blocks = std::min(out_count, available_blocks);
        size_t allocate_bytes = allocate_blocks * BLOCK_SIZE;
        char *next = cursor + allocate_bytes;

        if (arena.bump_cursor.compare_exchange_weak(cursor, next,
                                                    std::memory_order_relaxed,
                                                    std::memory_order_relaxed))
        {
            FreeBlock *head = reinterpret_cast<FreeBlock *>(cursor);
            FreeBlock *tail = head;

            for (size_t i = 1; i < allocate_blocks; ++i)
            {
                FreeBlock *block = reinterpret_cast<FreeBlock *>(cursor + i * BLOCK_SIZE);
                tail->next = block;
                tail = block;
            }

            tail->next = nullptr;
            *out_tail = tail;
            out_count = allocate_blocks;
            return head;
        }
    }
}

inline size_t ConcurrentMemoryPool::batch_allocate_from_arena(Arena &arena, FreeBlock **out_head, FreeBlock **out_tail, size_t count)
{
    FreeBlock *head = nullptr;
    FreeBlock *tail = nullptr;
    size_t fetched = 0;

    {
        std::lock_guard<std::mutex> lock(arena.mutex);

        if (arena.central_free_list)
        {
            head = arena.central_free_list;
            tail = head;
            fetched = 1;

            while (tail->next && fetched < count)
            {
                tail = tail->next;
                ++fetched;
            }

            // 更新中央链表
            arena.central_free_list = tail->next;
            tail->next = nullptr;
        }
    }

    if (fetched > 0)
    {
        *out_head = head;
        *out_tail = tail;
        return fetched;
    }

    size_t bump_count = count;
    head = batch_allocate_from_bump(arena, bump_count, &tail);
    if (!head)
    {
        return 0;
    }

    *out_head = head;
    *out_tail = tail;
    return bump_count;
}

inline void ConcurrentMemoryPool::batch_deallocate_to_arena(Arena &arena, FreeBlock *head, FreeBlock *tail, size_t /*count*/)
{
    if (!head || !tail)
        return;

    std::lock_guard<std::mutex> lock(arena.mutex);
    tail->next = arena.central_free_list;
    arena.central_free_list = head;
}

inline int ConcurrentMemoryPool::find_arena_index(void *ptr) const
{
    // 线性搜索（Arena 数量很少，最多 16 个）
    for (int i = 0; i < num_arenas_; ++i)
    {
        if (arenas_[i].contains(ptr))
        {
            return i;
        }
    }
    return -1;
}

inline int ConcurrentMemoryPool::get_thread_arena_index() const
{
#if HAS_LIBNUMA
    if (numa_available_)
    {
        int cpu = sched_getcpu();
        if (cpu >= 0)
        {
            int node = numa_node_of_cpu(cpu);
            if (node >= 0 && node < num_arenas_)
            {
                return node;
            }
        }
    }
#endif
    return 0;
}

inline ThreadLocalCache &ConcurrentMemoryPool::get_thread_local_cache()
{
    return tls_cache_;
}

inline void *ConcurrentMemoryPool::allocate_large(size_t bytes)
{
    if (bytes == 0)
    {
        return nullptr;
    }

    ThreadLocalCache &tls = tls_cache_;

    // 初始化 TLS（首次访问）
    if (!tls.pool) [[unlikely]]
    {
        tls.pool = this;
        tls.local_arena_index = get_thread_arena_index();
        tls.local_arena = &arenas_[tls.local_arena_index];
    }

    char *ptr = bump_allocate_from_arena(*tls.local_arena, bytes);
    if (ptr)
    {
        return ptr;
    }

    int start_index = (tls.local_arena_index + 1) % num_arenas_;
    for (int i = 0; i < num_arenas_; ++i)
    {
        int candidate = (start_index + i) % num_arenas_;
        if (candidate == tls.local_arena_index)
            continue;

        ptr = bump_allocate_from_arena(arenas_[candidate], bytes);
        if (ptr)
        {
            return ptr;
        }
    }

    std::cerr << "std::bad_alloc" << std::endl;
    std::exit(-1);
}

inline void *ConcurrentMemoryPool::allocate()
{
    ThreadLocalCache &tls = tls_cache_;

    // 初始化 TLS（首次访问）
    if (!tls.pool) [[unlikely]]
    {
        tls.pool = this;
        tls.local_arena_index = get_thread_arena_index();
        tls.local_arena = &arenas_[tls.local_arena_index];
    }

    // 1. 尝试从本地空闲栈分配
    if (tls.local_free_stack)
    {
        FreeBlock *block = tls.local_free_stack;
        tls.local_free_stack = block->next;
        tls.local_free_count--;
        return block;
    }

    // 2. 本地栈为空，从本地 Arena 批量获取
    FreeBlock *head, *tail;
    size_t fetched = batch_allocate_from_arena(*tls.local_arena, &head, &tail, BATCH_SIZE);
    if (fetched > 0)
    {
        // 第一个块立即返回
        FreeBlock *result = head;

        // 直接挂接剩余链表，避免额外遍历与反转
        tls.local_free_stack = head->next;
        result->next = nullptr;
        tls.local_free_count = fetched - 1;
        return result;
    }

    // 3. 本地 Arena 为空，尝试跨 Arena 窃取
    int start_index = (tls.local_arena_index + 1) % num_arenas_;
    for (int i = 0; i < num_arenas_; ++i)
    {
        int candidate = (start_index + i) % num_arenas_;
        if (candidate == tls.local_arena_index)
            continue;

        fetched = batch_allocate_from_arena(arenas_[candidate], &head, &tail, BATCH_SIZE);
        if (fetched > 0)
        {
            FreeBlock *result = head;

            // 与本地 refill 分支保持一致：直接挂接剩余链表
            tls.local_free_stack = head->next;
            result->next = nullptr;
            tls.local_free_count = fetched - 1;
            return result;
        }
    }

    // 4. 所有 Arena 都为空
    std::cerr << "std::bad_alloc" << std::endl;
    std::exit(-1);
}

inline void ConcurrentMemoryPool::deallocate(void *ptr)
{
    if (!ptr)
        return;

    ThreadLocalCache &tls = tls_cache_;

    // 初始化 TLS（首次访问）
    if (!tls.pool)
    {
        tls.pool = this;
        tls.local_arena_index = get_thread_arena_index();
        tls.local_arena = &arenas_[tls.local_arena_index];
    }

    FreeBlock *block = reinterpret_cast<FreeBlock *>(ptr);

    // 查找块属于哪个 Arena
    int arena_idx = find_arena_index(ptr);
    if (arena_idx < 0) [[unlikely]]
    {
        // 不属于任何 Arena，可能是无效指针
        return;
    }

    if (arena_idx == tls.local_arena_index)
    {
        // 本地 Arena
        if (tls.local_free_count < TLS_CAPACITY)
        {
            // 压入本地栈
            block->next = tls.local_free_stack;
            tls.local_free_stack = block;
            tls.local_free_count++;
        }
        else
        {
            // 栈已满，批量归还一半
            size_t half = TLS_CAPACITY / 2;
            FreeBlock *tail = tls.local_free_stack;
            for (size_t i = 1; i < half && tail->next; ++i)
            {
                tail = tail->next;
            }

            FreeBlock *remaining = tail->next;
            tail->next = nullptr;

            batch_deallocate_to_arena(*tls.local_arena, tls.local_free_stack, tail, half);

            tls.local_free_stack = remaining;
            tls.local_free_count = TLS_CAPACITY - half;

            // 将新块压入栈
            block->next = tls.local_free_stack;
            tls.local_free_stack = block;
            tls.local_free_count++;
        }
    }
    else
    {
        // 远程 Arena
        RemoteListHead &rl = tls.remote_lists[arena_idx];
        block->next = rl.head;
        rl.head = block;
        rl.count++;

        if (rl.count >= REMOTE_LIST_CAPACITY)
        {
            // 批量归还到远程 Arena
            FreeBlock *head = rl.head;
            // 找到尾节点
            FreeBlock *tail = head;
            while (tail->next)
            {
                tail = tail->next;
            }

            batch_deallocate_to_arena(arenas_[arena_idx], head, tail, rl.count);

            rl.head = nullptr;
            rl.count = 0;
        }
    }
}

inline int ConcurrentMemoryPool::get_current_arena_index() const
{
    return get_thread_arena_index();
}

inline void ConcurrentMemoryPool::reclaim_thread_cache()
{
    // TLS 的析构函数会自动处理
    // 这里可以显式触发清理
    ThreadLocalCache &tls = tls_cache_;
    if (tls.pool)
    {
        tls.~ThreadLocalCache();
        new (&tls) ThreadLocalCache();
    }
}

#endif // MEMORY_POOL_HPP
