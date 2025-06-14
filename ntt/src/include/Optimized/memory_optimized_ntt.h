#pragma once

#include "../config.h"
#include "../OpenMP_Barrett/ntt.h"
#include <immintrin.h>
#include <cstring>
#include <thread>
#include <iostream>

// 缓存行大小定义
#define CACHE_LINE_SIZE 64

// 内存对齐的数据结构
struct alignas(CACHE_LINE_SIZE) AlignedNTTData
{
    u64 data[CACHE_LINE_SIZE / sizeof(u64)];

    AlignedNTTData()
    {
        memset(data, 0, sizeof(data));
    }
};

// 线程局部存储，避免False Sharing
thread_local AlignedNTTData thread_local_buffer;

/**
 * 内存优化的NTT正变换
 * 主要优化：
 * 1. 缓存对齐的数据访问
 * 2. 避免False Sharing
 * 3. 优化内存访问模式
 * 4. 数据预取
 */
template <typename T>
inline void ntt_forward_memory_optimized(T *a, T n, T p, T omega)
{
    BarrettMod<T> mod(p, 100);

    // 预分配对齐的临时缓冲区
    size_t aligned_size = (n * sizeof(T) + CACHE_LINE_SIZE - 1) & ~(CACHE_LINE_SIZE - 1);
    T *aligned_a = (T *)_mm_malloc(aligned_size, CACHE_LINE_SIZE);

    // 复制数据到对齐缓冲区
    memcpy(aligned_a, a, n * sizeof(T));

    for (T mid = 1; mid < n; mid <<= 1)
    {
        T Wn = mod.pow(omega, (p - 1) / (mid << 1));

        // 计算每个线程处理的数据块大小，确保缓存对齐
        T block_size = (2 * mid + CACHE_LINE_SIZE / sizeof(T) - 1) & ~(CACHE_LINE_SIZE / sizeof(T) - 1);

#pragma omp parallel for schedule(static)
        for (T j = 0; j < n; j += (mid << 1))
        {
            // 数据预取优化
            if (j + (mid << 1) < n)
            {
                _mm_prefetch((const char *)&aligned_a[j + (mid << 1)], _MM_HINT_T0);
            }

            // 使用线程局部变量避免False Sharing
            T w = 1;

            // 内部循环展开优化
            T k = 0;
            for (; k + 3 < mid; k += 4, w = mod.mul(w, Wn))
            {
                // 批量处理4个元素，提高缓存利用率
                T x0 = aligned_a[j + k];
                T x1 = aligned_a[j + k + 1];
                T x2 = aligned_a[j + k + 2];
                T x3 = aligned_a[j + k + 3];

                T w1 = w;
                T w2 = mod.mul(w1, Wn);
                T w3 = mod.mul(w2, Wn);
                T w4 = mod.mul(w3, Wn);

                T y0 = mod.mul(w1, aligned_a[j + k + mid]);
                T y1 = mod.mul(w2, aligned_a[j + k + mid + 1]);
                T y2 = mod.mul(w3, aligned_a[j + k + mid + 2]);
                T y3 = mod.mul(w4, aligned_a[j + k + mid + 3]);

                aligned_a[j + k] = mod.add(x0, y0);
                aligned_a[j + k + 1] = mod.add(x1, y1);
                aligned_a[j + k + 2] = mod.add(x2, y2);
                aligned_a[j + k + 3] = mod.add(x3, y3);

                aligned_a[j + k + mid] = mod.sub(x0, y0);
                aligned_a[j + k + mid + 1] = mod.sub(x1, y1);
                aligned_a[j + k + mid + 2] = mod.sub(x2, y2);
                aligned_a[j + k + mid + 3] = mod.sub(x3, y3);

                w = w4;
            }

            // 处理剩余元素
            for (; k < mid; ++k, w = mod.mul(w, Wn))
            {
                T x = aligned_a[j + k];
                T y = mod.mul(w, aligned_a[j + k + mid]);
                aligned_a[j + k] = mod.add(x, y);
                aligned_a[j + k + mid] = mod.sub(x, y);
            }
        }
    }

    // 复制结果回原数组
    memcpy(a, aligned_a, n * sizeof(T));

    // 释放对齐的内存
    _mm_free(aligned_a);
}

/**
 * 内存优化的NTT逆变换
 */
template <typename T>
inline void ntt_backward_memory_optimized(T *a, T n, T p, T omega)
{
    BarrettMod<T> mod(p, 100);
    T omega_inv = mod.pow(omega, p - 2);

    // 使用相同的内存优化策略
    size_t aligned_size = (n * sizeof(T) + CACHE_LINE_SIZE - 1) & ~(CACHE_LINE_SIZE - 1);
    T *aligned_a = (T *)_mm_malloc(aligned_size, CACHE_LINE_SIZE);
    memcpy(aligned_a, a, n * sizeof(T));

    for (T mid = n >> 1; mid >= 1; mid >>= 1)
    {
        T Wn = mod.pow(omega_inv, (p - 1) / (mid << 1));

#pragma omp parallel for schedule(static)
        for (T j = 0; j < n; j += (mid << 1))
        {
            if (j + (mid << 1) < n)
            {
                _mm_prefetch((const char *)&aligned_a[j + (mid << 1)], _MM_HINT_T0);
            }

            T w = 1;
            for (T k = 0; k < mid; ++k, w = mod.mul(w, Wn))
            {
                T x = aligned_a[j + k];
                T y = aligned_a[j + k + mid];
                aligned_a[j + k] = mod.add(x, y);
                aligned_a[j + k + mid] = mod.mul(mod.sub(x, y), w);
            }
        }
    }

    // 归一化
    T n_inv = mod.pow(n, p - 2);
#pragma omp parallel for
    for (T i = 0; i < n; ++i)
    {
        aligned_a[i] = mod.mul(aligned_a[i], n_inv);
    }

    memcpy(a, aligned_a, n * sizeof(T));
    _mm_free(aligned_a);
}

/**
 * 内存优化的多项式乘法
 */
template <typename T>
inline void poly_multiply_ntt_memory_optimized(T *a, T *b, T *ab, T n, T p, T omega)
{
    T n_expanded = expand_n(2 * n - 1);

    // 使用缓存对齐的临时数组
    T *a_expanded = (T *)_mm_malloc(n_expanded * sizeof(T), CACHE_LINE_SIZE);
    T *b_expanded = (T *)_mm_malloc(n_expanded * sizeof(T), CACHE_LINE_SIZE);

    // 初始化并复制数据
    memset(a_expanded, 0, n_expanded * sizeof(T));
    memset(b_expanded, 0, n_expanded * sizeof(T));
    memcpy(a_expanded, a, n * sizeof(T));
    memcpy(b_expanded, b, n * sizeof(T));

    // 执行优化的NTT
    ntt_forward_memory_optimized(a_expanded, n_expanded, p, omega);
    ntt_forward_memory_optimized(b_expanded, n_expanded, p, omega);

    // 点乘优化
#pragma omp parallel for
    for (T i = 0; i < n_expanded; ++i)
    {
        a_expanded[i] = BarrettMod<T>(p).mul(a_expanded[i], b_expanded[i]);
    }

    // 逆变换
    ntt_backward_memory_optimized(a_expanded, n_expanded, p, omega);

    // 复制结果
    memcpy(ab, a_expanded, n_expanded * sizeof(T));

    // 清理内存
    _mm_free(a_expanded);
    _mm_free(b_expanded);
}

/**
 * 内存优化统计信息
 */
struct MemoryOptimizationStats
{
    size_t cache_aligned_allocations = 0;
    size_t prefetch_operations = 0;
    size_t loop_unroll_count = 0;

    void reset()
    {
        cache_aligned_allocations = 0;
        prefetch_operations = 0;
        loop_unroll_count = 0;
    }

    void print() const
    {
        std::cout << "[内存优化统计]" << std::endl;
        std::cout << "缓存对齐分配次数: " << cache_aligned_allocations << std::endl;
        std::cout << "数据预取操作次数: " << prefetch_operations << std::endl;
        std::cout << "循环展开次数: " << loop_unroll_count << std::endl;
    }
};

// 全局统计对象
extern MemoryOptimizationStats g_memory_stats;