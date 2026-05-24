#ifndef SCHEDULER_THREAD_POOL_HEADER
#define SCHEDULER_THREAD_POOL_HEADER

#include "definition.h"
#include "LayerQueues.h"
#include "MPMCRingQueue.h"
#include "NewKmerTree.h"

#include <memory>
#include <vector>
#include <thread>
#include <chrono>
#include <atomic>
#include <algorithm>
#include <utility>
#include <barrier>
#include <limits>

constexpr uint64_t TSL_HASH_TABLE_SIZE = 512 * 1024;

template <uint32_t N>
class SchedulerThreadPool final
{

    static constexpr uint32_t WORKER_QUEUE_CAPACITY = 4;
    static constexpr uint32_t FORMAL_WORKER_BATCH = 100;

    static constexpr uint32_t BACKLOG_PER_WORKER_TARGET = 2;
    static constexpr uint32_t MIN_BACKLOG_GATE_PRODUCING = 2;
    static constexpr uint32_t MIN_BACKLOG_GATE_DRAINING = 1;
    static constexpr uint32_t PRESSURE_WINDOW_BASE = 16;
    static constexpr uint32_t INVALID_DEPTH = MAX_DEPTH;

    // 以下是调度算法的参数
    static constexpr double SCORE_EMPTY_PENALTY = 0.20;
    static constexpr double SCORE_ACTIVE_DEPTH_BIAS = 0.03;
    static constexpr double SCORE_UPSTREAM_WEIGHT = 0.25;
    static constexpr double SCORE_BURST_BONUS = 1.2;
    static constexpr double SCORE_TREND_POS_WEIGHT = 1.1;
    static constexpr double SCORE_TREND_NEG_WEIGHT = 0.55;
    static constexpr double SCORE_SCALE_UP_THRESHOLD = 0.45;
    static constexpr double SCORE_SCALE_DOWN_THRESHOLD = -0.25;

    // 以下是worker线程状态机的参数
    static constexpr uint32_t SPIN_LIMIT_WAIT_DEPTH = 128;
    static constexpr uint32_t SPIN_LIMIT_BOUND_DEPTH = 64;
    static constexpr uint32_t MAX_BACKOFF_WAIT_DEPTH = 128;
    static constexpr uint32_t MAX_BACKOFF_BOUND_DEPTH = 64;
    static constexpr uint32_t RELEASE_ROUNDS_PRODUCING = 4;
    static constexpr uint32_t RELEASE_ROUNDS_DRAINING = 2;
    static constexpr uint32_t FORCE_LOCAL_STACK_EVERY_K_CANNOT_SWITCH = 64;
    static constexpr uint32_t DRAIN_EMPTY_CONFIRM_ROUNDS = 2;

    const uint32_t worker_count_;
    alignas(CACHE_LINE_SIZE) std::atomic<bool> running_{false};
    alignas(CACHE_LINE_SIZE) std::atomic<bool> stop_requested_{false};
    alignas(CACHE_LINE_SIZE) std::atomic<uint64_t> active_producer{0};

    KmerTree<N> *tree_ptr_ = nullptr;
    LayerQueues<N> *layer_queues_ptr_ = nullptr;
    std::thread scheduler_thread_;
    std::vector<std::unique_ptr<std::thread>> worker_threads_ptr_;
    std::vector<MPMCRingQueue<uint32_t, WORKER_QUEUE_CAPACITY>> worker_queues_;

    std::array<std::atomic<int>, MAX_DEPTH> real_active_workers;
    std::vector<std::atomic<uint32_t>> worker_last_depth;
    std::vector<std::atomic<uint8_t>> worker_is_bound;
    std::vector<std::atomic<uint32_t>> worker_current_depth;

// 测试模式
#ifdef TEST_MODE
    alignas(CACHE_LINE_SIZE) std::atomic<uint64_t> switch_success_count{0};
    alignas(CACHE_LINE_SIZE) std::atomic<uint64_t> force_local_stack_count{0};
    alignas(CACHE_LINE_SIZE) std::atomic<uint64_t> depth_release_count{0};
    alignas(CACHE_LINE_SIZE) std::atomic<uint64_t> wait_depth_yield_count{0};
    alignas(CACHE_LINE_SIZE) std::atomic<uint64_t> bound_depth_yield_count{0};
#endif

public:
    explicit SchedulerThreadPool(uint32_t worker_count, uint32_t producer_count, KmerTree<N> *tree_ptr, LayerQueues<N> *layer_queues_ptr)
        : worker_count_(worker_count), active_producer(producer_count), tree_ptr_(tree_ptr), layer_queues_ptr_(layer_queues_ptr), worker_queues_(worker_count - 1), worker_last_depth(worker_count - 1), worker_is_bound(worker_count - 1), worker_current_depth(worker_count - 1)
    {
        worker_threads_ptr_.reserve(worker_count_ - 1);
        for (uint32_t i = 0; i + 1 < worker_count_; i++)
        {
            worker_last_depth[i].store(0, std::memory_order_relaxed);
            worker_is_bound[i].store(0, std::memory_order_relaxed);
            worker_current_depth[i].store(INVALID_DEPTH, std::memory_order_relaxed);
        }
    }

    ~SchedulerThreadPool()
    {
#ifdef TEST_MODE
        std::cout << "Switch success count : " << switch_success_count.load(std::memory_order_relaxed) << std::endl;
        std::cout << "Force local stack count : " << force_local_stack_count.load(std::memory_order_relaxed) << std::endl;
        std::cout << "Depth release count : " << depth_release_count.load(std::memory_order_relaxed) << std::endl;
        std::cout << "Wait depth yield count : " << wait_depth_yield_count.load(std::memory_order_relaxed) << std::endl;
        std::cout << "Bound depth yield count : " << bound_depth_yield_count.load(std::memory_order_relaxed) << std::endl;
#endif
    }

    void start()
    {
        running_.store(true, std::memory_order_release);
        for (uint32_t i = 0; i < MAX_DEPTH; i++)
        {
            real_active_workers[i].store(0, std::memory_order_relaxed);
        }
        for (uint32_t i = 0; i + 1 < worker_count_; ++i)
        {
            worker_threads_ptr_.push_back(std::make_unique<std::thread>(&SchedulerThreadPool::worker_thread_loop, this, i));
        }
        scheduler_thread_ = std::thread(&SchedulerThreadPool::scheduler_thread_loop, this);
    }

    void join()
    {
        stop_requested_.store(true, std::memory_order_release);
        if (scheduler_thread_.joinable())
        {
            scheduler_thread_.join();
        }
        for (auto &t : worker_threads_ptr_)
        {
            if (t->joinable())
            {
                t->join();
            }
        }
    }

    void mark_producer_done()
    {
        active_producer.fetch_sub(1, std::memory_order_release);
    }

private:
    enum class WorkerState : uint8_t
    {
        WAIT_DEPTH = 0,
        BOUND_DEPTH = 1
    };

    struct WorkerLoopContext
    {
        WorkerState state = WorkerState::WAIT_DEPTH;
        uint32_t current_depth = 0;
        bool has_depth = false;
        uint32_t wait_spin_count = 0;
        uint32_t bound_spin_count = 0;
        uint32_t wait_backoff_iterations = 1;
        uint32_t bound_backoff_iterations = 1;
        uint32_t empty_rounds_at_depth = 0;
        uint32_t cannot_switch_count = 0;
    };

    bool try_bind_depth(const uint32_t worker_id, uint32_t &current_depth, bool &has_depth)
    {
        uint32_t next_depth = 0;
        if (!worker_queues_[worker_id].try_dequeue(next_depth))
        {
            return false;
        }

        if (next_depth >= MAX_DEPTH)
        {
            return false;
        }

        current_depth = next_depth;
        has_depth = true;
        worker_current_depth[worker_id].store(current_depth, std::memory_order_relaxed);
        real_active_workers[current_depth].fetch_add(1, std::memory_order_relaxed);
        worker_is_bound[worker_id].store(1, std::memory_order_release);
        return true;
    }

    bool try_switch_depth(const uint32_t worker_id, uint32_t &current_depth)
    {
        uint32_t next_depth = current_depth;
        if (!worker_queues_[worker_id].try_dequeue(next_depth))
        {
            return false;
        }

        if (next_depth >= MAX_DEPTH)
        {
            return false;
        }

        if (next_depth == current_depth)
        {
            return true;
        }

        real_active_workers[current_depth].fetch_sub(1, std::memory_order_relaxed);
        current_depth = next_depth;
        worker_current_depth[worker_id].store(current_depth, std::memory_order_relaxed);
        real_active_workers[current_depth].fetch_add(1, std::memory_order_relaxed);
#ifdef TEST_MODE
        switch_success_count.fetch_add(1, std::memory_order_relaxed);
#endif
        return true;
    }

    uint32_t process_batch_at_depth(const uint32_t depth, Task<N> &task)
    {
        auto queue = layer_queues_ptr_->get_queue(depth);
        uint32_t processed = 0;
        while (processed < FORMAL_WORKER_BATCH && queue->try_dequeue(task))
        {
            layer_queues_ptr_->decrease_size();
            tree_ptr_->thread_add_kmer(task);
            ++processed;
        }
        return processed;
    }

    bool maybe_force_process_local_stack(uint32_t &cannot_switch_count)
    {
        if (cannot_switch_count < FORCE_LOCAL_STACK_EVERY_K_CANNOT_SWITCH)
        {
            return false;
        }

        tree_ptr_->check_and_deal_with_local_stack();
        cannot_switch_count = 0;
#ifdef TEST_MODE
        force_local_stack_count.fetch_add(1, std::memory_order_relaxed);
#endif
        return true;
    }

    void release_depth_if_bound(const uint32_t worker_id, const uint32_t current_depth, bool &has_depth)
    {
        if (!has_depth)
        {
            return;
        }

        real_active_workers[current_depth].fetch_sub(1, std::memory_order_relaxed);
        has_depth = false;
        worker_current_depth[worker_id].store(INVALID_DEPTH, std::memory_order_relaxed);
        worker_is_bound[worker_id].store(0, std::memory_order_release);
#ifdef TEST_MODE
        depth_release_count.fetch_add(1, std::memory_order_relaxed);
#endif
    }

    void reset_bound_state_counters(WorkerLoopContext &context)
    {
        context.bound_spin_count = 0;
        context.bound_backoff_iterations = 1;
        context.empty_rounds_at_depth = 0;
    }

    bool are_all_depth_queues_empty() const
    {
        for (uint32_t depth = 0; depth < MAX_DEPTH; ++depth)
        {
            if (!layer_queues_ptr_->get_queue(depth)->empty())
            {
                return false;
            }
        }
        return true;
    }

    void apply_backoff(uint32_t &backoff_iterations, const uint32_t max_backoff, const bool is_wait_depth)
    {
        const uint32_t pause_count = std::min(backoff_iterations, max_backoff);
        for (uint32_t i = 0; i < pause_count; ++i)
        {
            cpu_relax();
        }

        if (backoff_iterations >= max_backoff)
        {
            std::this_thread::yield();
            backoff_iterations = 1;
#ifdef TEST_MODE
            if (is_wait_depth)
            {
                wait_depth_yield_count.fetch_add(1, std::memory_order_relaxed);
            }
            else
            {
                bound_depth_yield_count.fetch_add(1, std::memory_order_relaxed);
            }
#endif
            return;
        }

        backoff_iterations = std::min(backoff_iterations * 2, max_backoff);
    }

    void handle_wait_depth_state(const uint32_t worker_id, WorkerLoopContext &context)
    {
        if (try_bind_depth(worker_id, context.current_depth, context.has_depth))
        {
            context.state = WorkerState::BOUND_DEPTH;
            context.wait_spin_count = 0;
            context.wait_backoff_iterations = 1;
            reset_bound_state_counters(context);
            context.cannot_switch_count = 0;
            return;
        }

        ++context.wait_spin_count;
        if (context.wait_spin_count >= SPIN_LIMIT_WAIT_DEPTH)
        {
            context.wait_spin_count = 0;
            apply_backoff(context.wait_backoff_iterations, MAX_BACKOFF_WAIT_DEPTH, true);
        }
        else
        {
            for (int i = 0; i < 8; ++i)
            {
                cpu_relax();
            }
        }
    }

    void handle_bound_depth_state(const uint32_t worker_id, const bool draining_mode, Task<N> &task, WorkerLoopContext &context)
    {
        const uint32_t processed = process_batch_at_depth(context.current_depth, task);
        if (processed > 0)
        {
            reset_bound_state_counters(context);
        }

        if (try_switch_depth(worker_id, context.current_depth))
        {
            context.cannot_switch_count = 0;
            reset_bound_state_counters(context);
            return;
        }

        ++context.cannot_switch_count;
        if (maybe_force_process_local_stack(context.cannot_switch_count))
        {
            reset_bound_state_counters(context);
            return;
        }

        if (processed > 0)
        {
            return;
        }

        ++context.bound_spin_count;
        if (context.bound_spin_count < SPIN_LIMIT_BOUND_DEPTH)
        {
            cpu_relax();
            return;
        }

        context.bound_spin_count = 0;
        ++context.empty_rounds_at_depth;
        if (const uint32_t release_rounds = draining_mode ? RELEASE_ROUNDS_DRAINING : RELEASE_ROUNDS_PRODUCING;
            context.empty_rounds_at_depth >= release_rounds)
        {
            release_depth_if_bound(worker_id, context.current_depth, context.has_depth);
            context.state = WorkerState::WAIT_DEPTH;
            context.empty_rounds_at_depth = 0;
            context.cannot_switch_count = 0;
            return;
        }

        apply_backoff(context.bound_backoff_iterations, MAX_BACKOFF_BOUND_DEPTH, false);
    }

    void drain_all_depths_and_local_stack(Task<N> &task)
    {
        uint32_t stable_empty_rounds = 0;
        while (stable_empty_rounds < DRAIN_EMPTY_CONFIRM_ROUNDS)
        {
            bool found_queue_work = false;
            for (int32_t depth = static_cast<int32_t>(MAX_DEPTH - 1); depth >= 0; --depth)
            {
                auto queue = layer_queues_ptr_->get_queue(static_cast<uint32_t>(depth));
                while (queue->try_dequeue(task))
                {
                    layer_queues_ptr_->decrease_size();
                    tree_ptr_->thread_add_kmer(task);
                    found_queue_work = true;
                }
            }

            tree_ptr_->check_and_deal_with_local_stack();

            if (!found_queue_work && are_all_depth_queues_empty())
            {
                ++stable_empty_rounds;
            }
            else
            {
                stable_empty_rounds = 0;
            }
        }
    }

    void worker_thread_loop(const uint32_t worker_id)
    {
        WorkerLoopContext context;
        bool draining_mode = active_producer.load(std::memory_order_acquire) == 0;
        Task<N> task;

        while (!stop_requested_.load(std::memory_order_relaxed) || layer_queues_ptr_->size() > 0)
        {
            if (!draining_mode && active_producer.load(std::memory_order_acquire) == 0)
            {
                draining_mode = true;
            }

            if (context.state == WorkerState::WAIT_DEPTH)
            {
                handle_wait_depth_state(worker_id, context);
                continue;
            }

            handle_bound_depth_state(worker_id, draining_mode, task, context);
        }

        release_depth_if_bound(worker_id, context.current_depth, context.has_depth);
        drain_all_depths_and_local_stack(task);

#ifdef TEST_MODE
        SpinLock::flush_spin_loops_for_current_thread();
#endif
    }

    double get_score(uint32_t depth, const std::array<uint32_t, MAX_DEPTH> &queue_size,
                     const std::array<double, MAX_DEPTH> &pressure, double pressure_diff)
    {
        double score = pressure[depth];

        if (queue_size[depth] == 0)
        {
            score -= SCORE_EMPTY_PENALTY; // 空队惩罚
        }
        else
        {
            score += depth * SCORE_ACTIVE_DEPTH_BIAS; // 按照深度有加分
        }

        if (depth > 0 && queue_size[depth - 1] > 0)
        {
            score += pressure[depth - 1] * SCORE_UPSTREAM_WEIGHT; // 受上层压力影响，因为很快回到下层
        }

        if (pressure[depth] > 0.9)
        {
            score += SCORE_BURST_BONUS; // 突发流量加分
        }

        if (pressure_diff > 0.0)
        {
            score += pressure_diff * SCORE_TREND_POS_WEIGHT; // 压力变化趋势加分，响应压力上升
        }
        else
        {
            score += pressure_diff * SCORE_TREND_NEG_WEIGHT; // 压力变化趋势减分，响应压力下降
        }

        return score;
    }

    void scheduler_thread_loop()
    {
        const uint32_t available_workers = worker_count_ - 1;
        const uint32_t max_quota = std::max<uint32_t>(1, available_workers * 3 / 4);

        std::array<uint32_t, MAX_DEPTH> min_threads_limit;
        std::array<uint32_t, MAX_DEPTH> soft_max_limit;
        std::array<uint32_t, MAX_DEPTH> hard_max_limit;
        std::array<uint32_t, MAX_DEPTH> max_threads_limit;
        std::array<uint32_t, MAX_DEPTH> dynamic_max_threads;

        for (uint32_t i = 0; i < MAX_DEPTH; i++)
        {
            min_threads_limit[i] = 1;
            if (i == 0)
            {
                soft_max_limit[i] = std::min<uint32_t>(
                    std::max<uint32_t>(2, available_workers / 4s),
                    max_quota);
                hard_max_limit[i] = max_quota;
            }
            else
            {
                soft_max_limit[i] = std::min<uint32_t>(2 * (1 << i), max_quota);
                hard_max_limit[i] = max_quota;
            }
            max_threads_limit[i] = soft_max_limit[i];
            dynamic_max_threads[i] = min_threads_limit[i];
        }

        std::array<std::pair<double, uint32_t>, MAX_DEPTH> scores;
        std::array<double, MAX_DEPTH> prev_pressure = {0.0};

        // 核心监控与调度循环
        while (!stop_requested_.load(std::memory_order_relaxed) || layer_queues_ptr_->size() > 0)
        {
            std::array<uint32_t, MAX_DEPTH> queue_size;
            std::array<double, MAX_DEPTH> fill_ratio;
            std::array<double, MAX_DEPTH> pressure;
            std::array<double, MAX_DEPTH> pressure_diff;
            const uint32_t pressure_window = std::max<uint32_t>(PRESSURE_WINDOW_BASE, available_workers * BACKLOG_PER_WORKER_TARGET);
            for (uint32_t i = 0; i < MAX_DEPTH; i++)
            {
                auto queue = layer_queues_ptr_->get_queue(i);
                const uint32_t q = static_cast<uint32_t>(queue->size());
                queue_size[i] = q;
                fill_ratio[i] = static_cast<double>(q) / static_cast<double>(TASK_QUEUE_CAPACITY);
                pressure[i] = std::min(1.0, static_cast<double>(q) / static_cast<double>(pressure_window));
                pressure_diff[i] = pressure[i] - prev_pressure[i];
                prev_pressure[i] = pressure[i];
            }

            for (uint32_t i = 0; i < MAX_DEPTH; i++)
            {
                const bool should_expand_cap =
                    fill_ratio[i] > 0.85 ||
                    (i == 0 && (pressure[i] > 0.85 || pressure_diff[i] > 0.12));
                const bool should_shrink_cap =
                    fill_ratio[i] < 0.6 &&
                    (i != 0 || (pressure[i] < 0.55 && pressure_diff[i] <= 0.0));

                // 1. 如果水压严重爆表，天花板临时打开，允许向硬顶靠拢
                if (should_expand_cap && max_threads_limit[i] < hard_max_limit[i])
                {
                    max_threads_limit[i]++;
                }
                // 2. 如果水压恢复正常（甚至低于平时），天花板慢慢降回软顶
                else if (should_shrink_cap && max_threads_limit[i] > soft_max_limit[i])
                {
                    max_threads_limit[i]--;
                }

                // 超限恢复校验：确保动态分配人数被收口
                if (dynamic_max_threads[i] > max_threads_limit[i])
                {
                    dynamic_max_threads[i] = max_threads_limit[i];
                }

                double score = get_score(i, queue_size, pressure, pressure_diff[i]);

                scores[i] = {score, i};

                // === 构造非线性扩容步长函数 ===
                int32_t step = 0;

                // 1. 基于绝对压力(score)决定基础步调
                if (score > SCORE_SCALE_UP_THRESHOLD)
                {
                    step = 1 + static_cast<int32_t>((score - SCORE_SCALE_UP_THRESHOLD) * 2.0);
                }
                else if (score < SCORE_SCALE_DOWN_THRESHOLD)
                {
                    step = -1 - static_cast<int32_t>((SCORE_SCALE_DOWN_THRESHOLD - score) * 1.5);
                }

                // 2. 叠加瞬时水流冲击率 (高能预警) 赋予额外补偿步长
                if (pressure_diff[i] > 0.08)
                {
                    step += static_cast<int32_t>(pressure_diff[i] * 10.0);
                }
                else if (pressure_diff[i] < -0.10 && step < 0)
                {
                    step -= 1;
                }

                // === 安全地施加该步长 ===
                if (step > 0)
                {
                    uint32_t increase = static_cast<uint32_t>(step);
                    dynamic_max_threads[i] = std::min<uint32_t>(dynamic_max_threads[i] + increase, max_threads_limit[i]);
                }
                else if (step < 0)
                {
                    uint32_t shrink = static_cast<uint32_t>(-step);
                    if (dynamic_max_threads[i] > min_threads_limit[i] + shrink)
                    {
                        dynamic_max_threads[i] -= shrink;
                    }
                    else
                    {
                        dynamic_max_threads[i] = min_threads_limit[i];
                    }
                }
            }

            std::sort(scores.begin(), scores.end(), [](const auto &a, const auto &b)
                      { return a.first > b.first; });

            const bool producers_done = active_producer.load(std::memory_order_acquire) == 0;
            const uint32_t min_backlog_gate = producers_done ? MIN_BACKLOG_GATE_DRAINING : MIN_BACKLOG_GATE_PRODUCING;

            std::array<uint32_t, MAX_DEPTH> target_threads;
            for (uint32_t depth = 0; depth < MAX_DEPTH; ++depth)
            {
                const uint32_t backlog = queue_size[depth];
                if (backlog < min_backlog_gate)
                {
                    target_threads[depth] = 0;
                    continue;
                }

                uint32_t target = (backlog + BACKLOG_PER_WORKER_TARGET - 1) / BACKLOG_PER_WORKER_TARGET;
                target = std::max<uint32_t>(target, min_threads_limit[depth]);
                target = std::min<uint32_t>(target, dynamic_max_threads[depth]);
                target_threads[depth] = target;
            }

            std::array<int32_t, MAX_DEPTH> planned_active;
            for (uint32_t depth = 0; depth < MAX_DEPTH; ++depth)
            {
                planned_active[depth] = real_active_workers[depth].load(std::memory_order_relaxed);
            }

            for (uint32_t worker_id = 0; worker_id < available_workers; ++worker_id)
            {
                if (worker_queues_[worker_id].empty())
                {
                    continue;
                }

                const uint32_t pending_depth = worker_last_depth[worker_id].load(std::memory_order_relaxed);
                if (pending_depth >= MAX_DEPTH)
                {
                    continue;
                }

                if (worker_is_bound[worker_id].load(std::memory_order_acquire) != 0)
                {
                    const uint32_t current_depth = worker_current_depth[worker_id].load(std::memory_order_relaxed);
                    if (current_depth < MAX_DEPTH && current_depth != pending_depth)
                    {
                        if (planned_active[current_depth] > 0)
                        {
                            planned_active[current_depth]--;
                        }
                        planned_active[pending_depth]++;
                    }
                }
                else
                {
                    planned_active[pending_depth]++;
                }
            }

            auto can_assign_depth = [&](const uint32_t depth) -> bool
            {
                if (depth >= MAX_DEPTH)
                {
                    return false;
                }

                if (queue_size[depth] < min_backlog_gate)
                {
                    return false;
                }

                return planned_active[depth] < static_cast<int32_t>(target_threads[depth]);
            };

            auto pick_best_deficit_depth = [&]() -> uint32_t
            {
                for (const auto &s : scores)
                {
                    const uint32_t depth = s.second;
                    if (can_assign_depth(depth))
                    {
                        return depth;
                    }
                }
                return INVALID_DEPTH;
            };

            for (uint32_t worker_id = 0; worker_id < available_workers; worker_id++)
            {
                if (worker_is_bound[worker_id].load(std::memory_order_acquire) != 0)
                {
                    continue;
                }

                if (!worker_queues_[worker_id].empty())
                {
                    continue;
                }

                const uint32_t last_depth = worker_last_depth[worker_id].load(std::memory_order_relaxed);
                uint32_t selected_depth = INVALID_DEPTH;

                if (can_assign_depth(last_depth))
                {
                    selected_depth = last_depth;
                }
                else
                {
                    selected_depth = pick_best_deficit_depth();
                }

                if (selected_depth != INVALID_DEPTH &&
                    worker_is_bound[worker_id].load(std::memory_order_acquire) == 0 &&
                    worker_queues_[worker_id].empty() &&
                    worker_queues_[worker_id].try_enqueue(selected_depth)) [[likely]]
                {
                    worker_last_depth[worker_id].store(selected_depth, std::memory_order_relaxed);
                    planned_active[selected_depth]++;
                }
            }

            auto is_burst_depth = [&](const uint32_t depth) -> bool
            {
                return fill_ratio[depth] > 0.75 ||
                       pressure[depth] > 0.85 ||
                       pressure_diff[depth] > 0.12;
            };

            auto should_migrate_depth = [&](const uint32_t depth) -> bool
            {
                if (depth >= MAX_DEPTH || target_threads[depth] == 0)
                {
                    return false;
                }

                const int32_t deficit = static_cast<int32_t>(target_threads[depth]) - planned_active[depth];
                if (deficit <= 0)
                {
                    return false;
                }

                return is_burst_depth(depth) || planned_active[depth] == 0 || deficit >= 2;
            };

            auto can_donate_depth = [&](const uint32_t depth) -> bool
            {
                if (depth >= MAX_DEPTH)
                {
                    return false;
                }

                if (pressure[depth] > 0.65)
                {
                    return false;
                }

                return planned_active[depth] > static_cast<int32_t>(target_threads[depth]);
            };

            auto find_donor_worker = [&](const uint32_t dst) -> int32_t
            {
                int32_t best_worker = -1;
                double best_pressure = 2.0;
                int32_t best_excess = std::numeric_limits<int32_t>::min();

                for (uint32_t worker_id = 0; worker_id < available_workers; ++worker_id)
                {
                    if (worker_is_bound[worker_id].load(std::memory_order_acquire) == 0)
                    {
                        continue;
                    }

                    if (!worker_queues_[worker_id].empty())
                    {
                        continue;
                    }

                    const uint32_t src = worker_current_depth[worker_id].load(std::memory_order_relaxed);
                    if (src >= MAX_DEPTH || src == dst || !can_donate_depth(src))
                    {
                        continue;
                    }

                    const int32_t excess = planned_active[src] - static_cast<int32_t>(target_threads[src]);
                    if (pressure[src] < best_pressure || (pressure[src] == best_pressure && excess > best_excess))
                    {
                        best_pressure = pressure[src];
                        best_excess = excess;
                        best_worker = static_cast<int32_t>(worker_id);
                    }
                }

                return best_worker;
            };

            uint32_t migration_budget = std::max<uint32_t>(1, available_workers / 8);
            for (const auto &s : scores)
            {
                if (migration_budget == 0)
                {
                    break;
                }

                const uint32_t dst = s.second;
                if (!should_migrate_depth(dst))
                {
                    continue;
                }

                while (migration_budget > 0 && can_assign_depth(dst))
                {
                    const int32_t donor = find_donor_worker(dst);
                    if (donor < 0)
                    {
                        break;
                    }

                    const uint32_t donor_id = static_cast<uint32_t>(donor);
                    if (worker_is_bound[donor_id].load(std::memory_order_acquire) == 0 || !worker_queues_[donor_id].empty())
                    {
                        break;
                    }

                    const uint32_t src = worker_current_depth[donor_id].load(std::memory_order_relaxed);
                    if (src >= MAX_DEPTH || src == dst || !can_donate_depth(src))
                    {
                        break;
                    }

                    if (worker_queues_[donor_id].try_enqueue(dst)) [[likely]]
                    {
                        worker_last_depth[donor_id].store(dst, std::memory_order_relaxed);
                        planned_active[src]--;
                        planned_active[dst]++;
                        --migration_budget;
                    }
                    else
                    {
                        break;
                    }
                }
            }

            std::this_thread::sleep_for(std::chrono::nanoseconds(500));
        }
    }
};

#endif // SCHEDULING_THREAD_POOL_HEADER
