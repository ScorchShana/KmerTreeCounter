#ifndef MPMC_RING_QUEUE_HEADER
#define MPMC_RING_QUEUE_HEADER

#include "definition.h"

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <optional>
#include <thread>
#include <type_traits>
#include <utility>



template <typename T, std::size_t Capacity>
class MPMCRingQueue
{

	enum State : std::uint8_t
	{
		EMPTY = 0,
		STORING = 1,
		STORED = 2,
		LOADING = 3
	};

public:
	static_assert((Capacity & (Capacity - 1)) == 0, "Capacity must be power of 2");
	static_assert(Capacity >= CACHE_LINE_SIZE / sizeof(std::atomic<State>), "Capacity must be at least the number of states per cache line");
	static_assert(std::is_trivially_copyable_v<T>, "T must be trivially copyable");

	explicit MPMCRingQueue() noexcept:
		datas_(Capacity), states_(Capacity)
 	{
		for (std::size_t i = 0; i < Capacity; ++i)
		{
			states_[i].store(EMPTY, std::memory_order_relaxed);
		}
	}
	

	inline bool enqueue(const T &item) noexcept
	{
		std::size_t pos = 0;
		
		constexpr uint32_t SPIN_LIMIT = 256;
		uint32_t spin_count = 1;
		while (!try_reserve_tail(pos))
		{
			for (uint32_t i = 0; i < spin_count; i++)
			{
				cpu_relax();
			}
			spin_count = std::min(spin_count * 2, SPIN_LIMIT);
			if (spin_count >= SPIN_LIMIT)
			{
				std::this_thread::yield();
				spin_count = 1;
			}
		}

		std::atomic<State> &state = states_[get_state_index(pos)];
		wait_and_claim_for_store(state);

		std::memcpy(static_cast<void *>(std::addressof(datas_[get_data_index(pos)])),
					static_cast<const void *>(std::addressof(item)),
					sizeof(T));
		// Publish data: memcpy above happens-before consumers that acquire STORED.
		state.store(STORED, std::memory_order_release);
		return true;
	}

	inline bool try_enqueue(const T &item) noexcept
	{
		std::size_t pos = 0;
		if (!try_reserve_tail(pos))
		{
			return false;
		}

		std::atomic<State> &state = states_[get_state_index(pos)];
		wait_and_claim_for_store(state);

		std::memcpy(static_cast<void *>(std::addressof(datas_[get_data_index(pos)])),
					static_cast<const void *>(std::addressof(item)),
					sizeof(T));
		// Publish data: memcpy above happens-before consumers that acquire STORED.
		state.store(STORED, std::memory_order_release);
		return true;
	}

	inline bool dequeue(T &out) noexcept
	{
		constexpr uint32_t SPIN_LIMIT = 256;
		uint32_t spin_count = 1;
		while(!try_dequeue(out))
		{
			for(uint32_t i = 0; i < spin_count; i++)
			{
				cpu_relax();
			}
			spin_count = std::min(spin_count * 2, SPIN_LIMIT);
			if (spin_count >= SPIN_LIMIT)
			{
				std::this_thread::yield();
				spin_count = 1;
			}
		}
		return true;
	}

	inline bool try_dequeue(T &out) noexcept
	{
		std::size_t pos = 0;
		if (!try_reserve_head(pos))
		{
			return false;
		}

		std::atomic<State> &state = states_[get_state_index(pos)];
		wait_and_claim_for_load(state);
		out = datas_[get_data_index(pos)];
		// Recycle slot: read above happens-before the next producer that acquires EMPTY.
		state.store(EMPTY, std::memory_order_release);
		return true;
	}

	inline bool try_pop(T &out) noexcept
	{
		return try_dequeue(out);
	}

	inline bool empty() const noexcept
	{
		const std::size_t head = head_.load(std::memory_order_relaxed);
		const std::size_t tail = tail_.load(std::memory_order_relaxed);
		return head == tail;
	}

	inline bool full() const noexcept
	{
		const std::size_t head = head_.load(std::memory_order_relaxed);
		const std::size_t tail = tail_.load(std::memory_order_relaxed);
		return (tail - head) >= Capacity;
	}

	inline std::size_t size() const noexcept
	{
		const std::size_t head = head_.load(std::memory_order_relaxed);
		const std::size_t tail = tail_.load(std::memory_order_relaxed);
		return tail - head;
	}

private:
	static constexpr std::size_t kMask = Capacity - 1;

	// 生产者侧：尝试在全局 tail 上预留一个写入位置（返回到 pos）。
	// 成功后表示当前线程独占了逻辑下标 pos 对应的槽位写权限（后续还需等待该槽位 state=EMPTY）。
	inline bool try_reserve_tail(std::size_t &pos) noexcept
	{
		// 从当前 tail 快照开始竞争。
		pos = tail_.load(std::memory_order_relaxed);
		for (;;)
		{
			// 条件 (pos - head) >= Capacity 表示环形队列“逻辑已满”：
			// [head, pos) 已经占满 Capacity 个元素，没有可预留的新写位置。
			if (const std::size_t head = head_.load(std::memory_order_acquire); (pos - head) >= Capacity)
			{
				return false;
			}
			// CAS 成功：把 tail 从 pos 推进到 pos+1，当前线程拿到写票据 pos。
			// CAS 失败：pos 会被更新为最新 tail，继续重试。
			if (tail_.compare_exchange_weak(pos, pos + 1, std::memory_order_relaxed, std::memory_order_relaxed))
			{
				return true;
			}
		}
	}

	// 消费者侧：尝试在全局 head 上预留一个读取位置（返回到 pos）。
	// 成功后表示当前线程独占了逻辑下标 pos 对应的槽位读权限（后续还需等待该槽位 state=STORED）。
	inline bool try_reserve_head(std::size_t &pos) noexcept
	{
		// 从当前 head 快照开始竞争。
		pos = head_.load(std::memory_order_relaxed);
		for (;;)
		{
			// 条件 (tail - pos) == 0 表示队列“逻辑为空”：
			// 没有任何尚未被消费的元素可读。
			if (const std::size_t tail = tail_.load(std::memory_order_acquire); (tail - pos) == 0)
			{
				return false;
			}
			// CAS 成功：把 head 从 pos 推进到 pos+1，当前线程拿到读票据 pos。
			// CAS 失败：pos 会被更新为最新 head，继续重试。
			if (head_.compare_exchange_weak(pos, pos + 1, std::memory_order_relaxed, std::memory_order_relaxed))
			{
				return true;
			}
		}
	}

	// 针对“某个已预留写位置”的槽位状态机：
	// 仅当槽位为 EMPTY 时，生产者才能把它原子地占用为 STORING，随后写入 datas_ 并发布为 STORED。
	inline static void wait_and_claim_for_store(std::atomic<State> &state) noexcept
	{
		State expected = EMPTY;
		// CAS 条件：expected==EMPTY 时成功，表示槽位可写；
		// 失败则说明槽位仍在被前一轮使用（可能是 STORING/STORED/LOADING），继续自旋等待。
		while (!state.compare_exchange_weak(expected, STORING, std::memory_order_acquire, std::memory_order_relaxed))
		{
			expected = EMPTY;
			cpu_relax();
		}
	}

	// 针对“某个已预留读位置”的槽位状态机：
	// 仅当槽位为 STORED 时，消费者才能把它原子地占用为 LOADING，随后读取 datas_ 并回收为 EMPTY。
	inline static void wait_and_claim_for_load(std::atomic<State> &state) noexcept
	{
		State expected = STORED;
		// CAS 条件：expected==STORED 时成功，表示该槽位已有完整数据可读；
		// 失败则说明数据尚未发布或被其他线程处理中，继续自旋等待。
		while (!state.compare_exchange_weak(expected, LOADING, std::memory_order_acquire, std::memory_order_relaxed))
		{
			expected = STORED;
			cpu_relax();
		}
	}

	inline static std::size_t get_state_index(std::size_t pos) noexcept
	{
		const std::size_t pos_in_cycle = pos & kMask;
		constexpr std::size_t state_size = sizeof(std::atomic<State>);
		constexpr std::size_t states_per_cache_line = CACHE_LINE_SIZE / state_size;
		constexpr std::size_t states_per_column = Capacity / states_per_cache_line;
		return (pos_in_cycle % states_per_column) * states_per_cache_line + (pos_in_cycle / states_per_column);
	}

	inline static std::size_t get_data_index(std::size_t pos) noexcept
	{
		return pos & kMask;
	}

	alignas(CACHE_LINE_SIZE) std::vector<T> datas_;
	alignas(CACHE_LINE_SIZE) std::vector<std::atomic<State>> states_;
	alignas(CACHE_LINE_SIZE) std::atomic<std::size_t> head_{0};
	alignas(CACHE_LINE_SIZE) std::atomic<std::size_t> tail_{0};
};

#endif // MPMC_RING_QUEUE_HEADER
