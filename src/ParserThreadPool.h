#ifndef PARSER_THREAD_POOL_HEADER
#define PARSER_THREAD_POOL_HEADER

#include "definition.h"
#include "MPMCRingQueue.h"
#include "NewKmerTree.h"
#include "ConcurrentMemoryPool.h"
#include "RingMemoryPool.h"
#include "FastqParser.h"
#include "SchedulerThreadPool.h"

#include <array>
#include <atomic>
#include <cstdint>
#include <memory>
#include <thread>
#include <vector>
#include <utility>
#include <chrono>

template <uint32_t N>
class ParserThreadPool final
{
public:
    explicit ParserThreadPool(const int in_k, KmerTree<N> *tree_ptr,
                              ConcurrentMemoryPool *pool_ptr,
                              RingMemoryPool<RING_MEMORY_POOL_CAPACITY> *ring_pool_ptr,
                              SchedulerThreadPool<N> *task_thread_pool_ptr,
                              uint32_t worker_count)
        : k(in_k),
          tree_(tree_ptr),
          pool_(pool_ptr),
          ring_pool_(ring_pool_ptr),
          task_thread_pool_(task_thread_pool_ptr),
          worker_count_(worker_count == 0 ? 1u : worker_count)
    {
        threads_ptr_.reserve(worker_count_);
    }

    void start()
    {
        for (uint32_t i = 0; i < worker_count_; ++i)
        {
            threads_ptr_.push_back(std::make_unique<std::thread>([&]
                                                                 {
				FastqParser<N> parser(k, ring_pool_, tree_, pool_);
				parser.parse_and_push();
				//pool_->flush_all_thread_caches_to_central();
				pool_->reclaim_thread_cache();
				total_read_kmer += parser.get_total_read_kmer();
                SpinLock::flush_spin_loops_for_current_thread();
                task_thread_pool_->mark_producer_done(); }));
        }
    }

    void join()
    {
        for (auto &t : threads_ptr_)
        {
            if (t->joinable())
            {
                t->join();
            }
        }
    }

    inline uint64_t get_total_read_kmer() const
    {
        return total_read_kmer.load(std::memory_order_acquire);
    }

    inline void add_total_read_kmer(uint64_t count) noexcept
    {
        total_read_kmer.fetch_add(count, std::memory_order_acq_rel);
    }

private:
    int k;
    KmerTree<N> *tree_ = nullptr;
    ConcurrentMemoryPool *pool_ = nullptr;
    RingMemoryPool<RING_MEMORY_POOL_CAPACITY> *ring_pool_ = nullptr;
    std::vector<std::unique_ptr<std::thread>> threads_ptr_;
    SchedulerThreadPool<N> *task_thread_pool_ = nullptr;
    const uint32_t worker_count_;
    std::atomic<uint64_t> total_read_kmer = 0;
};

#endif // PARSER_THREAD_POOL_HEADER