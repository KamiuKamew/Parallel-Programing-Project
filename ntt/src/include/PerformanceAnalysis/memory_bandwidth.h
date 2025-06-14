#pragma once

#include <vector>
#include <thread>
#include <atomic>
#include <cstring>
#include "../config.h"

/**
 * @brief 内存带宽测试结果结构体
 */
struct MemoryBandwidthResult
{
    double bandwidth_gbps;   // 内存带宽 (GB/s)
    double latency_ns;       // 平均访问延迟 (纳秒)
    double time_ms;          // 总执行时间 (毫秒)
    size_t bytes_accessed;   // 访问的字节数
    int num_threads;         // 线程数
    double efficiency_ratio; // 相对于单线程的效率比例
};

/**
 * @brief 内存带宽竞争分析结果
 */
struct MemoryContentionAnalysis
{
    std::vector<MemoryBandwidthResult> sequential_results; // 顺序访问结果
    std::vector<MemoryBandwidthResult> random_results;     // 随机访问结果
    std::vector<MemoryBandwidthResult> stride_results;     // 跨步访问结果
    double peak_bandwidth_gbps;                            // 峰值带宽
    int optimal_thread_count;                              // 最优线程数
};

/**
 * @brief 顺序内存访问测试
 */
void sequential_memory_access(volatile uint64_t *data, size_t size, size_t iterations)
{
    for (size_t iter = 0; iter < iterations; ++iter)
    {
        for (size_t i = 0; i < size; i += 8)
        { // 64字节cache line，每次跳8个uint64_t
            volatile uint64_t temp = data[i];
            data[i] = temp + 1; // 读写操作
        }
    }
}

/**
 * @brief 随机内存访问测试
 */
void random_memory_access(volatile uint64_t *data, size_t size, size_t iterations)
{
    // 预生成随机索引避免测试期间的随机数生成开销
    std::vector<size_t> indices(iterations);
    for (size_t i = 0; i < iterations; ++i)
    {
        indices[i] = (i * 1103515245 + 12345) % size; // 简单的LCG
    }

    for (size_t iter = 0; iter < iterations; ++iter)
    {
        size_t idx = indices[iter];
        volatile uint64_t temp = data[idx];
        data[idx] = temp + 1;
    }
}

/**
 * @brief 跨步内存访问测试
 */
void stride_memory_access(volatile uint64_t *data, size_t size, size_t stride, size_t iterations)
{
    for (size_t iter = 0; iter < iterations; ++iter)
    {
        for (size_t i = 0; i < size; i += stride)
        {
            volatile uint64_t temp = data[i];
            data[i] = temp + 1;
        }
    }
}

/**
 * @brief 测试特定线程数下的内存带宽
 */
MemoryBandwidthResult measure_memory_bandwidth(int num_threads, size_t memory_size_mb = 256,
                                               size_t iterations = 1000, int access_pattern = 0)
{
    MemoryBandwidthResult result;
    result.num_threads = num_threads;

    size_t total_size = memory_size_mb * 1024 * 1024 / sizeof(uint64_t); // 转换为uint64_t数量
    std::vector<uint64_t> memory_pool(total_size, 0);
    volatile uint64_t *data = memory_pool.data();

    // 计算每个线程的数据范围
    size_t chunk_size = total_size / num_threads;
    std::vector<std::thread> threads;
    std::atomic<bool> start_flag{false};

    TIMER_START();

    // 创建线程
    for (int t = 0; t < num_threads; ++t)
    {
        size_t start_idx = t * chunk_size;
        size_t end_idx = (t == num_threads - 1) ? total_size : (t + 1) * chunk_size;
        size_t thread_size = end_idx - start_idx;

        threads.emplace_back([&start_flag, data, start_idx, thread_size, iterations, access_pattern]()
                             {
            // 等待同步启动
            while (!start_flag.load()) {
                std::this_thread::yield();
            }
            
            volatile uint64_t* thread_data = data + start_idx;
            
            switch (access_pattern) {
                case 0: // 顺序访问
                    sequential_memory_access(thread_data, thread_size, iterations);
                    break;
                case 1: // 随机访问
                    random_memory_access(thread_data, thread_size, iterations);
                    break;
                case 2: // 跨步访问
                    stride_memory_access(thread_data, thread_size, 64, iterations);  // 64-stride
                    break;
            } });
    }

    // 同步启动所有线程
    start_flag.store(true);

    // 等待所有线程完成
    for (auto &thread : threads)
    {
        thread.join();
    }

    TIMER_END();
    result.time_ms = TIMER_ELAPSED();

    // 计算带宽统计
    result.bytes_accessed = total_size * sizeof(uint64_t) * iterations * 2; // 读写各一次
    result.bandwidth_gbps = (result.bytes_accessed / 1e9) / (result.time_ms / 1000.0);
    result.latency_ns = (result.time_ms * 1e6) / (total_size * iterations); // 平均每次访问的延迟

    return result;
}

/**
 * @brief 综合内存竞争分析
 */
MemoryContentionAnalysis analyze_memory_contention(const std::vector<int> &thread_counts,
                                                   size_t memory_size_mb = 256,
                                                   size_t iterations = 1000)
{
    MemoryContentionAnalysis analysis;
    analysis.peak_bandwidth_gbps = 0.0;
    analysis.optimal_thread_count = 1;

    // 测试不同访问模式
    for (size_t pattern = 0; pattern < 3; ++pattern)
    {
        std::vector<MemoryBandwidthResult> *results_ptr;

        switch (pattern)
        {
        case 0:
            results_ptr = &analysis.sequential_results;
            break;
        case 1:
            results_ptr = &analysis.random_results;
            break;
        case 2:
            results_ptr = &analysis.stride_results;
            break;
        }

        // 获得单线程基准
        MemoryBandwidthResult baseline = measure_memory_bandwidth(1, memory_size_mb, iterations, pattern);
        baseline.efficiency_ratio = 1.0;
        results_ptr->push_back(baseline);

        // 测试不同线程数
        for (int threads : thread_counts)
        {
            if (threads == 1)
                continue; // 已经测试过

            MemoryBandwidthResult result = measure_memory_bandwidth(threads, memory_size_mb, iterations, pattern);
            result.efficiency_ratio = result.bandwidth_gbps / (baseline.bandwidth_gbps * threads);
            results_ptr->push_back(result);

            // 更新峰值带宽和最优线程数
            if (result.bandwidth_gbps > analysis.peak_bandwidth_gbps)
            {
                analysis.peak_bandwidth_gbps = result.bandwidth_gbps;
                analysis.optimal_thread_count = threads;
            }
        }
    }

    return analysis;
}

/**
 * @brief Cache友好性测试 - 测试false sharing影响
 */
struct CacheFriendlinessResult
{
    double aligned_time_ms;        // 缓存对齐版本时间
    double unaligned_time_ms;      // 未对齐版本时间
    double false_sharing_overhead; // false sharing开销比例
    int num_threads;
};

CacheFriendlinessResult test_false_sharing_impact(int num_threads, size_t iterations = 10000000)
{
    CacheFriendlinessResult result;
    result.num_threads = num_threads;

    // 测试对齐版本 (每个线程一个cache line)
    {
        struct alignas(64) AlignedCounter
        {
            volatile uint64_t value = 0;
            char padding[64 - sizeof(uint64_t)];
        };

        std::vector<AlignedCounter> aligned_counters(num_threads);
        std::vector<std::thread> threads;

        TIMER_START();
        for (int t = 0; t < num_threads; ++t)
        {
            threads.emplace_back([&aligned_counters, t, iterations]()
                                 {
                for (size_t i = 0; i < iterations; ++i) {
                    aligned_counters[t].value++;
                } });
        }

        for (auto &thread : threads)
        {
            thread.join();
        }
        TIMER_END();
        result.aligned_time_ms = TIMER_ELAPSED();
    }

    // 测试未对齐版本 (可能发生false sharing)
    {
        std::vector<uint64_t> unaligned_counters(num_threads, 0);
        std::vector<std::thread> threads;

        TIMER_START();
        for (int t = 0; t < num_threads; ++t)
        {
            threads.emplace_back([&unaligned_counters, t, iterations]()
                                 {
                for (size_t i = 0; i < iterations; ++i) {
                    unaligned_counters[t]++;
                } });
        }

        for (auto &thread : threads)
        {
            thread.join();
        }
        TIMER_END();
        result.unaligned_time_ms = TIMER_ELAPSED();
    }

    result.false_sharing_overhead = result.unaligned_time_ms / result.aligned_time_ms;

    return result;
}