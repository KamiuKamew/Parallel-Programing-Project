#pragma once

#include "unified.h"
#include "thread_pool_unified.h"
#include "../general/utils.h"
#include <thread>
#include <algorithm>
#include <memory>

/**
 * @brief 统一的Pthread并行NTT前向变换
 *
 * @param a 多项式系数数组（已在计算域中）
 * @param n 多项式长度
 * @param p 模数
 * @param omega_domain 原根（计算域）
 * @param w_domain 初始旋转因子（计算域）
 * @param mod_impl 模运算实现
 */
inline void ntt_forward_unified_pthread(u32 *a, u32 n, u32 p, u32 omega_domain, u32 w_domain, ModBase *mod_impl)
{
    size_t num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0)
        num_threads = 1;

    UnifiedThreadPool pool(num_threads);

    for (u32 mid = 1; mid < n; mid <<= 1)
    {
        u32 Wn_domain = mod_impl->pow(omega_domain, (p - 1) / (mid << 1));
        u32 total_blocks = n / (mid << 1);

        for (size_t t = 0; t < num_threads; ++t)
        {
            u32 start_block = total_blocks * t / num_threads;
            u32 end_block = total_blocks * (t + 1) / num_threads;

            if (start_block < end_block)
            {
                pool.enqueue([=, &mod_impl]()
                             {
                    for (u32 b = start_block; b < end_block; ++b)
                    {
                        u32 j = b * (mid << 1);
                        u32 w = w_domain;
                        for (u32 k = 0; k < mid; ++k, w = mod_impl->mul(w, Wn_domain))
                        {
                            u32 x = a[j + k];
                            u32 y = mod_impl->mul(w, a[j + k + mid]);
                            a[j + k] = mod_impl->add(x, y);
                            a[j + k + mid] = mod_impl->sub(x, y);
                        }
                    } });
            }
        }

        pool.wait();
    }
}

/**
 * @brief 统一的Pthread并行NTT逆变换
 *
 * @param a 频域系数数组（已在计算域中）
 * @param n 多项式长度
 * @param p 模数
 * @param omega_domain 原根（计算域）
 * @param w_domain 初始旋转因子（计算域）
 * @param mod_impl 模运算实现
 */
inline void ntt_inverse_unified_pthread(u32 *a, u32 n, u32 p, u32 inv_omega_domain, u32 w_domain, ModBase *mod_impl)
{
    size_t num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0)
        num_threads = 1;

    UnifiedThreadPool pool(num_threads);

    for (u32 mid = n >> 1; mid > 0; mid >>= 1)
    {
        u32 Wn_domain = mod_impl->pow(inv_omega_domain, (p - 1) / (mid << 1));
        u32 total_blocks = n / (mid << 1);

        for (size_t t = 0; t < num_threads; ++t)
        {
            u32 start_block = total_blocks * t / num_threads;
            u32 end_block = total_blocks * (t + 1) / num_threads;

            if (start_block < end_block)
            {
                pool.enqueue([=, &mod_impl]()
                             {
                    for (u32 b = start_block; b < end_block; ++b)
                    {
                        u32 j = b * (mid << 1);
                        u32 w = w_domain;
                        for (u32 k = 0; k < mid; ++k, w = mod_impl->mul(w, Wn_domain))
                        {
                            u32 x = a[j + k];
                            u32 y = a[j + k + mid];
                            a[j + k] = mod_impl->add(x, y);
                            a[j + k + mid] = mod_impl->mul(w, mod_impl->sub(x, y));
                        }
                    } });
            }
        }

        pool.wait();
    }

    // 归一化
    u32 inv_n = mod_impl->inv(mod_impl->to_compute_domain(n));
    for (u32 i = 0; i < n; ++i)
        a[i] = mod_impl->mul(a[i], inv_n);
}

/**
 * @brief 统一的Pthread并行多项式乘法
 *
 * @param a 第一个多项式系数
 * @param b 第二个多项式系数
 * @param ab 输出多项式系数
 * @param n 多项式长度
 * @param p 模数
 * @param omega 原根
 * @param mod_impl 模运算实现
 */
inline void poly_multiply_unified_pthread(u32 *a, u32 *b, u32 *ab, u32 n, u32 p, u32 omega, ModBase *mod_impl)
{
    u32 n_expanded = n << 1;
    u32 omega_domain = mod_impl->to_compute_domain(omega);
    u32 w_domain = mod_impl->to_compute_domain(1);
    u32 inv_omega_domain = mod_impl->inv(omega_domain);

    // 扩展数组
    u32 *a_expanded = expand_a(a, n, n_expanded);
    u32 *b_expanded = expand_a(b, n, n_expanded);

    // 先比特翻转
    bit_reverse_permute(a_expanded, n_expanded);
    bit_reverse_permute(b_expanded, n_expanded);

    // 再转换到计算域
    mod_impl->array_to_compute_domain(a_expanded, n_expanded);
    mod_impl->array_to_compute_domain(b_expanded, n_expanded);

    // 前向NTT
    ntt_forward_unified_pthread(a_expanded, n_expanded, p, omega_domain, w_domain, mod_impl);
    ntt_forward_unified_pthread(b_expanded, n_expanded, p, omega_domain, w_domain, mod_impl);

    // 点乘
    for (u32 i = 0; i < n_expanded; ++i)
        a_expanded[i] = mod_impl->mul(a_expanded[i], b_expanded[i]);

    // 逆向NTT
    ntt_inverse_unified_pthread(a_expanded, n_expanded, p, inv_omega_domain, w_domain, mod_impl);

    // 转换回普通域
    mod_impl->array_from_compute_domain(a_expanded, n_expanded);
    // 输出前再做一次比特翻转
    bit_reverse_permute(a_expanded, n_expanded);
    for (u32 i = 0; i < n; ++i)
        ab[i] = a_expanded[i];

    delete[] a_expanded;
    delete[] b_expanded;
}
