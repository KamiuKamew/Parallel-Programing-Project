#pragma once

#include <omp.h>
#include <vector>
#include <chrono>
#include "../config.h"
#include "../Barrett/op.h"

/**
 * @brief 同步开销测试结果结构体
 */
struct SyncOverheadResult
{
    double total_time_ms;       // 总执行时间(毫秒)
    double computation_time_ms; // 纯计算时间(毫秒)
    double sync_time_ms;        // 同步等待时间(毫秒)
    double sync_overhead_ratio; // 同步开销占比 (sync_time / total_time)
    int num_threads;            // 线程数
    int num_barriers;           // barrier次数
    double avg_barrier_time_ms; // 平均单次barrier时间(毫秒)
};

/**
 * @brief 带详细同步时间统计的NTT前向变换
 */
template <typename T>
SyncOverheadResult ntt_forward_with_sync_analysis(T *a, T n, T p, T omega, int num_threads)
{
    SyncOverheadResult result;
    result.num_threads = num_threads;
    result.num_barriers = 0;
    result.sync_time_ms = 0.0;

    BarrettMod<T> mod(p, 100);
    omp_set_num_threads(num_threads);

    TIMER_START();

    for (T mid = 1; mid < n; mid <<= 1)
    {
        T Wn = mod.pow(omega, (p - 1) / (mid << 1));

        // 记录计算开始时间
        double computation_start = omp_get_wtime();

#pragma omp parallel
        {
            // 每个线程记录自己的同步时间
            double thread_sync_time = 0.0;

#pragma omp for
            for (T j = 0; j < n; j += (mid << 1))
            {
                T w = 1;
                for (T k = 0; k < mid; ++k, w = mod.mul(w, Wn))
                {
                    T x = a[j + k];
                    T y = mod.mul(w, a[j + k + mid]);
                    a[j + k] = mod.add(x, y);
                    a[j + k + mid] = mod.sub(x, y);
                }
            }

            // 测量 barrier 时间
            double barrier_start = omp_get_wtime();
#pragma omp barrier
            double barrier_end = omp_get_wtime();
            thread_sync_time += (barrier_end - barrier_start);

#pragma omp atomic
            result.sync_time_ms += thread_sync_time * 1000.0;
        }

        result.num_barriers++;
    }

    TIMER_END();
    result.total_time_ms = TIMER_ELAPSED();

    // 由于每个线程都累加了同步时间，需要除以线程数得到平均值
    result.sync_time_ms /= num_threads;
    result.computation_time_ms = result.total_time_ms - result.sync_time_ms;
    result.sync_overhead_ratio = result.sync_time_ms / result.total_time_ms;
    result.avg_barrier_time_ms = result.sync_time_ms / result.num_barriers;

    return result;
}

/**
 * @brief 测试不同线程数下的同步开销
 */
template <typename T>
std::vector<SyncOverheadResult> analyze_sync_overhead_scaling(T n, T p, T omega,
                                                              const std::vector<int> &thread_counts,
                                                              int iterations = 5)
{
    std::vector<SyncOverheadResult> results;

    for (int num_threads : thread_counts)
    {
        SyncOverheadResult avg_result;
        avg_result.num_threads = num_threads;
        avg_result.total_time_ms = 0.0;
        avg_result.computation_time_ms = 0.0;
        avg_result.sync_time_ms = 0.0;
        avg_result.sync_overhead_ratio = 0.0;
        avg_result.avg_barrier_time_ms = 0.0;
        avg_result.num_barriers = 0;

        for (int iter = 0; iter < iterations; ++iter)
        {
            // 准备测试数据
            std::vector<T> test_data(n);
            for (T i = 0; i < n; ++i)
            {
                test_data[i] = i % p;
            }

            // 执行测试
            SyncOverheadResult iter_result = ntt_forward_with_sync_analysis(
                test_data.data(), n, p, omega, num_threads);

            // 累加结果
            avg_result.total_time_ms += iter_result.total_time_ms;
            avg_result.computation_time_ms += iter_result.computation_time_ms;
            avg_result.sync_time_ms += iter_result.sync_time_ms;
            avg_result.sync_overhead_ratio += iter_result.sync_overhead_ratio;
            avg_result.avg_barrier_time_ms += iter_result.avg_barrier_time_ms;
            avg_result.num_barriers = iter_result.num_barriers; // 所有迭代相同
        }

        // 计算平均值
        avg_result.total_time_ms /= iterations;
        avg_result.computation_time_ms /= iterations;
        avg_result.sync_time_ms /= iterations;
        avg_result.sync_overhead_ratio /= iterations;
        avg_result.avg_barrier_time_ms /= iterations;

        results.push_back(avg_result);
    }

    return results;
}

/**
 * @brief 测试静态vs动态调度的同步开销差异
 */
template <typename T>
struct SchedulingOverheadResult
{
    SyncOverheadResult static_result;
    SyncOverheadResult dynamic_result;
    SyncOverheadResult guided_result;
    double dynamic_overhead_ratio; // dynamic相对于static的额外开销
    double guided_overhead_ratio;  // guided相对于static的额外开销
};

template <typename T>
SchedulingOverheadResult<T> compare_scheduling_overhead(T n, T p, T omega, int num_threads)
{
    SchedulingOverheadResult<T> result;

    BarrettMod<T> mod(p, 100);
    omp_set_num_threads(num_threads);

    // 测试静态调度
    {
        std::vector<T> test_data(n);
        for (T i = 0; i < n; ++i)
            test_data[i] = i % p;

        TIMER_START();
        for (T mid = 1; mid < n; mid <<= 1)
        {
            T Wn = mod.pow(omega, (p - 1) / (mid << 1));
#pragma omp parallel for schedule(static)
            for (T j = 0; j < n; j += (mid << 1))
            {
                T w = 1;
                for (T k = 0; k < mid; ++k, w = mod.mul(w, Wn))
                {
                    T x = test_data[j + k];
                    T y = mod.mul(w, test_data[j + k + mid]);
                    test_data[j + k] = mod.add(x, y);
                    test_data[j + k + mid] = mod.sub(x, y);
                }
            }
        }
        TIMER_END();
        result.static_result.total_time_ms = TIMER_ELAPSED();
        result.static_result.num_threads = num_threads;
    }

    // 测试动态调度
    {
        std::vector<T> test_data(n);
        for (T i = 0; i < n; ++i)
            test_data[i] = i % p;

        TIMER_START();
        for (T mid = 1; mid < n; mid <<= 1)
        {
            T Wn = mod.pow(omega, (p - 1) / (mid << 1));
#pragma omp parallel for schedule(dynamic, 1)
            for (T j = 0; j < n; j += (mid << 1))
            {
                T w = 1;
                for (T k = 0; k < mid; ++k, w = mod.mul(w, Wn))
                {
                    T x = test_data[j + k];
                    T y = mod.mul(w, test_data[j + k + mid]);
                    test_data[j + k] = mod.add(x, y);
                    test_data[j + k + mid] = mod.sub(x, y);
                }
            }
        }
        TIMER_END();
        result.dynamic_result.total_time_ms = TIMER_ELAPSED();
        result.dynamic_result.num_threads = num_threads;
    }

    // 测试guided调度
    {
        std::vector<T> test_data(n);
        for (T i = 0; i < n; ++i)
            test_data[i] = i % p;

        TIMER_START();
        for (T mid = 1; mid < n; mid <<= 1)
        {
            T Wn = mod.pow(omega, (p - 1) / (mid << 1));
#pragma omp parallel for schedule(guided)
            for (T j = 0; j < n; j += (mid << 1))
            {
                T w = 1;
                for (T k = 0; k < mid; ++k, w = mod.mul(w, Wn))
                {
                    T x = test_data[j + k];
                    T y = mod.mul(w, test_data[j + k + mid]);
                    test_data[j + k] = mod.add(x, y);
                    test_data[j + k + mid] = mod.sub(x, y);
                }
            }
        }
        TIMER_END();
        result.guided_result.total_time_ms = TIMER_ELAPSED();
        result.guided_result.num_threads = num_threads;
    }

    // 计算额外开销比例
    result.dynamic_overhead_ratio = result.dynamic_result.total_time_ms / result.static_result.total_time_ms;
    result.guided_overhead_ratio = result.guided_result.total_time_ms / result.static_result.total_time_ms;

    return result;
}