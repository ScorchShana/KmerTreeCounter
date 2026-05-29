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
class ParserThreadPool
{

    int k;
    ConcurrentMemoryPool *pool = nullptr;
    SPMCRingMemoryPool<READER_PARSER_RING_MEMORY_POOL_CAPACITY> *reader_parser_ring_pool;
    RingMemoryPool<PARSER_CLASSIFIER_RING_MEMORY_POOL_CAPACITY> *parser_classifier_ring_pool;
    std::vector<std::unique_ptr<std::thread>> threads_ptr;
    const uint32_t parser_count;
    std::atomic<uint64_t> total_read_kmer = 0;

public:
#ifdef TEST_MODE
    std::atomic<uint64_t> producer_enqueue_spin_time{0};
    std::atomic<uint64_t> producer_dequeue_spin_time{0};
    std::atomic<uint64_t> consumer_enqueue_spin_time{0};
    std::atomic<uint64_t> consumer_dequeue_spin_time{0};
#endif

    explicit ParserThreadPool(const int in_k,
                              ConcurrentMemoryPool *pool_ptr,
                              SPMCRingMemoryPool<READER_PARSER_RING_MEMORY_POOL_CAPACITY> *in_reader_parser_ring_pool_ptr,
                              RingMemoryPool<PARSER_CLASSIFIER_RING_MEMORY_POOL_CAPACITY> *in_parser_classifier_ring_pool_ptr,
                              uint32_t in_parser_count)
        : k(in_k),
          pool(pool_ptr),
          reader_parser_ring_pool(in_reader_parser_ring_pool_ptr),
          parser_classifier_ring_pool(in_parser_classifier_ring_pool_ptr),
          parser_count(in_parser_count)
    {
        threads_ptr.reserve(parser_count);
    }

    void start()
    {
        for (uint32_t i = 0; i < parser_count; ++i)
        {
            threads_ptr.push_back(std::make_unique<std::thread>([&]
                                                                {
                                                                    std::this_thread::sleep_for(std::chrono::nanoseconds(20));
                                                                    FastqParser<N> parser(k, reader_parser_ring_pool, parser_classifier_ring_pool);
                                                                    parser.parse_and_push();
                                                                    total_read_kmer += parser.get_total_read_kmer();
                                                                    parser_classifier_ring_pool->producer_set_finished();

#ifdef TEST_MODE
                                                                    producer_enqueue_spin_time += parser.producer_enqueue_spin_time;
                                                                    producer_dequeue_spin_time += parser.producer_dequeue_spin_time;
                                                                    consumer_enqueue_spin_time += parser.consumer_enqueue_spin_time;
                                                                    consumer_dequeue_spin_time += parser.consumer_dequeue_spin_time;
#endif
                                                                }));
        }
    }

    void join()
    {
        for (auto &t : threads_ptr)
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
};

#endif // PARSER_THREAD_POOL_HEADER