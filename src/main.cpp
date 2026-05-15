#include "NewKmerTree.h"
#include "SchedulerThreadPool.h"
#include "ConcurrentMemoryPool.h"
#include "MPMCRingQueue.h"
#include "GetKmer.h"
#include "FastqReader.h"
#include "ParserThreadPool.h"
#include "ExportWriter.h"

#include <chrono>
#include <iostream>
#include <cstdint>
#include <string>
#include <memory>
#include <limits>
#include <format>


uint32_t k_len; // k-mer 的长度 (例如: 31, 41)
uint32_t n_thread; // 允许使用的总线程数
uint64_t memory_limit; // 全局最大内存限制，单位为 GB

int main(int argc, char *argv[])
{

    if (argc < 5 || argc > 10)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <fastq_file> <k_len> <n_thread> <memory_limit_gb> [map_capacity] [min_count] [max_count] [parser_threads]" << std::endl;
        return 1;
    }
    bool enable_prefault = false;
    uint32_t parser_threads_override = 0;
    min_count = 1;
    max_count = std::numeric_limits<uint32_t>::max();

    std::string filename;
    try
    {
        filename = argv[1];
        k_len = std::stoul(argv[2]);
        n_thread = std::stoul(argv[3]);
        memory_limit = std::stoull(argv[4]);

        if (argc >= 6)
        {
            kmer_concurrent_hash_map_capacity = std::stoul(argv[5]);
        }
        if (argc >= 7)
        {
            min_count = std::stoul(argv[6]);
        }
        if (argc >= 8)
        {
            max_count = std::stoul(argv[7]);
        }
        if (argc >= 9)
        {
            parser_threads_override = std::stoul(argv[8]);
            if (parser_threads_override <= 0 || parser_threads_override >= n_thread)
            {
                std::cerr << "parser_threads must be > 0 and < n_thread" << std::endl;
                return 1;
            }
        }

        if (n_thread < 3)
        {
            throw std::invalid_argument("n_thread must be >= 3");
        }

        std::cout << "Input parameters:" << std::endl;
        std::cout << "  Fastq file: " << filename << std::endl;
        std::cout << "  k-mer length: " << k_len << std::endl;
        std::cout << "  Thread count: " << n_thread << std::endl;
        std::cout << "  Memory limit (GB): " << memory_limit << std::endl;
        std::cout << "  Map capacity: " << kmer_concurrent_hash_map_capacity << std::endl;
        std::cout << "  Min count: " << min_count << std::endl;
        std::cout << "  Max count: " << max_count << std::endl;
        if (parser_threads_override > 0)
        {
            std::cout << "  Parser threads override: " << parser_threads_override << std::endl;
        }
    }
    catch (const std::exception &)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <fastq_file> <k_len> <n_thread> <memory_limit_gb> [map_capacity] [min_count] [max_count] [parser_threads]" << std::endl;
        return 1;
    }

    if (kmer_concurrent_hash_map_capacity <= 1 || kmer_concurrent_hash_map_capacity >= 16ULL * 1024 * 1024 || max_count < min_count)
    {
        std::cerr << "Usage: " << argv[0]
                  << " <fastq_file> <k_len> <n_thread> <memory_limit_gb> [map_capacity] [min_count] [max_count] [parser_threads]" << std::endl;
        return 1;
    }

    constexpr uint32_t N1 = 2;

    // 根据预算计算分给 Worker 的解析(Parser)线程和任务(Tasker)线程的数量
    // Worker 预算去除了主线程(Reader)和导出线程(ExportWriter) 
    const uint32_t worker_budget = n_thread - 2;
    uint32_t parser_num = 1 + (worker_budget - 1) / (1 + TASK_PARSER_RATIO);
    if (parser_threads_override > 0)
    {
        parser_num = parser_threads_override;
    }
    if (parser_num == 0)
    {
        parser_num = 1;
    }
    if (parser_num >= worker_budget)
    {
        parser_num = worker_budget;
    }
    uint32_t tasker_num = worker_budget - parser_num;
    tasker_num = tasker_num > 0 ? tasker_num : 1;

    std::cout << "Thread split:" << std::endl;
    std::cout << "  parser threads: " << parser_num << std::endl;
    std::cout << "  task threads: " << tasker_num << std::endl;

    const auto init_start = std::chrono::steady_clock::now();

    // 初始化层级队列，用于在树的不同深度间传递任务
    auto layer_queues = std::make_shared<LayerQueues<N1>>();
    // 初始化解析器环形内存池，管理 Parser 解析后的 k-mer 数据块
    auto ring_pool = std::make_shared<RingMemoryPool<RING_MEMORY_POOL_CAPACITY>>(RING_MEMORY_POOL_BLOCK_SIZE, 1, parser_num);
    // 初始化导出用的环形内存池，管理低频 k-mer 的导出数据块
    auto export_ring_pool = std::make_shared<RingMemoryPool<EXPORT_RING_MEMORY_POOL_CAPACITY>>(EXPORT_RING_MEMORY_POOL_BLOCK_SIZE, 1, 1);
    // 初始化全局并发内存池，用于节点分配、哈希表等
    auto pool = std::make_shared<ConcurrentMemoryPool>(memory_limit * 1024ULL * 1024ULL * 1024ULL, n_thread);

    // 初始化 k-mer 字典树(KmerTree)的核心结构，整合前述多个组件
    auto tree = std::make_shared<KmerTree<N1>>(k_len, pool.get(), layer_queues.get(), export_ring_pool.get());

    // 初始化并构建 Tasker 线程池，负责消费层级队列并将 k-mer 路由到深层节点 / 哈希表
    auto task_thread_pool = std::make_shared<SchedulerThreadPool<N1>>(tasker_num, parser_num, tree.get(), layer_queues.get());
    // 初始化并构建 Parser 线程池，负责消费 FASTQ 读取器产生的数据，提取 k-mer 进行初步布隆过滤
    auto parser_thread_pool = std::make_shared<ParserThreadPool<N1>>(k_len, tree.get(), pool.get(), ring_pool.get(), task_thread_pool.get(), parser_num);

    // 初始化导出写入器，用于将低频 k-mer 单线程安全落盘
    auto export_writer = std::make_shared<ExportWriter<N1>>(export_ring_pool.get());

    // 初始化最终排干(drain)阶段所需的初始任务
    layer_queues->initialize_final_drain_queue(tree->root_nodes);

    // 初始化 FASTQ 读取器，将大文件分块读取并送入 ring_pool 用作流水线起点
    FastqReader<N1> reader(filename, k_len, FASTQ_FILE_CHUNK_SIZE, ring_pool.get());

    const auto init_end = std::chrono::steady_clock::now();

    const auto mid_start = std::chrono::steady_clock::now();

    // 启动多线程流水线的各个工作线程
    parser_thread_pool->start();
    task_thread_pool->start();
    export_writer->start();

    // 在主线程中运行读取器，读取文件块并推送给 Parser
    reader.read();

    const auto read_end = std::chrono::steady_clock::now();

    // 阻塞等待 Parser 线程将所有文件块解析并生成 k-mer 完成
    parser_thread_pool->join();
    const auto parser_end = std::chrono::steady_clock::now();
    // 阻塞等待 Tasker 线程消费完所有的分发任务
    task_thread_pool->join();

    const auto task_end = std::chrono::steady_clock::now();

    const auto mid_end = std::chrono::steady_clock::now();

    const auto final_start = std::chrono::steady_clock::now();

    // Final drain 阶段：多线程并行遍历整个字典树，将在节点中暂存但未下发的 k-mers 全部合并到全局哈希表中
    tree->final_drain_parallel(n_thread - 1);

    // 标记所有的缓存已经发送完毕，向导出环形队列发送 nullptr 或退出标志
    tree->mark_finish_export();

    // 阻塞等待导出器将所有低频 k-mer 都安全写入磁盘完成
    export_writer->join();

    std::cout << "Total read k-mer count: " << parser_thread_pool->get_total_read_kmer() << std::endl;

    const auto final_end = std::chrono::steady_clock::now();

    const auto init_elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(init_end - init_start).count();
    const auto mid_elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(mid_end - mid_start).count();
    const auto final_elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(final_end - final_start).count();
    const auto total_elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(final_end - init_start).count();
    const auto parse_task_elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(task_end - parser_end).count();
    const auto read_parse_elapsed_us = std::chrono::duration_cast<std::chrono::microseconds>(parser_end - read_end).count();

    std::cout << "Init elapsed us: " << init_elapsed_us << std::endl;
    std::cout << "Mid elapsed us: " << mid_elapsed_us << std::endl;
    std::cout << "Between Read and Parse elapsed us: " << read_parse_elapsed_us << std::endl;
    std::cout << "Between Parse and Task elapsed us: " << parse_task_elapsed_us << std::endl;
    std::cout << "Final elapsed us: " << final_elapsed_us << std::endl;
    std::cout << "Total elapsed us: " << total_elapsed_us << std::endl;

#ifdef TEST_MODE
    SpinLock::flush_spin_loops_for_current_thread();
    std::cout << "SpinLock spin_loops: " << SpinLock::spin_loops() << std::endl;

    /*for (uint64_t i = 0; i < (1ULL << (2 * ROOT_BASES)); ++i)
    {
        std::cout << "Root prefix " << std::format("{:0{}d}", i, 3) << " count: " << root_counts[i].load(std::memory_order_relaxed) << std::endl;
    }*/
#endif

    return 0;
}