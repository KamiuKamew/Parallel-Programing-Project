#pragma once

#include "unified.h"
#include "../general/utils.h"
#include <omp.h>

/**
 * @brief 统一的OpenMP并行NTT前向变换
 *
 * @param a 多项式系数数组（已在计算域中）
 * @param n 多项式长度
 * @param p 模数
 * @param omega_domain 原根（计算域）
 * @param w_domain 初始旋转因子（计算域）
 * @param mod_impl 模运算实现
 */
inline void ntt_forward_unified_omp(u32 *a, u32 n, u32 p, u32 omega_domain, u32 w_domain, ModBase *mod_impl)
{
    for (u32 mid = 1; mid < n; mid <<= 1)
    {
        u32 Wn_domain = mod_impl->pow(omega_domain, (p - 1) / (mid << 1));

// OpenMP并行化j循环
#pragma omp parallel for
        for (u32 j = 0; j < n; j += (mid << 1))
        {
            u32 _w_domain = w_domain;
            for (u32 k = 0; k < mid; ++k)
            {
                u32 x_domain = a[j + k];
                u32 y_domain = mod_impl->mul(_w_domain, a[j + k + mid]);
                a[j + k] = mod_impl->add(x_domain, y_domain);
                a[j + k + mid] = mod_impl->sub(x_domain, y_domain);
                _w_domain = mod_impl->mul(_w_domain, Wn_domain);
            }
        }
    }
}

/**
 * @brief 统一的OpenMP并行NTT逆变换
 *
 * @param a 频域系数数组（已在计算域中）
 * @param n 多项式长度
 * @param p 模数
 * @param omega_domain 原根（计算域）
 * @param w_domain 初始旋转因子（计算域）
 * @param mod_impl 模运算实现
 */
inline void ntt_inverse_unified_omp(u32 *a, u32 n, u32 p, u32 omega_domain, u32 w_domain, ModBase *mod_impl)
{
    u32 omega_inv_domain = mod_impl->inv(omega_domain);

    for (u32 mid = n >> 1; mid > 0; mid >>= 1)
    {
        u32 Wn_domain = mod_impl->pow(omega_inv_domain, (p - 1) / (mid << 1));

// OpenMP并行化j循环
#pragma omp parallel for
        for (u32 j = 0; j < n; j += (mid << 1))
        {
            u32 _w_domain = w_domain;
            for (u32 k = 0; k < mid; ++k)
            {
                u32 x_domain = a[j + k];
                u32 y_domain = a[j + k + mid];
                a[j + k] = mod_impl->add(x_domain, y_domain);
                a[j + k + mid] = mod_impl->mul(_w_domain, mod_impl->sub(x_domain, y_domain));
                _w_domain = mod_impl->mul(_w_domain, Wn_domain);
            }
        }
    }

    // 归一化：乘以n的逆元
    u32 inv_n_domain = mod_impl->inv(mod_impl->to_compute_domain(n));
#pragma omp parallel for
    for (u32 i = 0; i < n; ++i)
        a[i] = mod_impl->mul(a[i], inv_n_domain);
}

/**
 * @brief 统一的OpenMP并行多项式乘法
 *
 * @param a 第一个多项式的系数
 * @param b 第二个多项式的系数
 * @param ab 结果多项式的系数
 * @param n 多项式长度
 * @param p 模数
 * @param omega 原根
 * @param mod_impl 模运算实现
 */
inline void poly_multiply_unified_omp(u32 *a, u32 *b, u32 *ab, u32 n, u32 p, u32 omega, ModBase *mod_impl)
{
    // 扩展到2*n-1的下一个2的幂
    u32 n_expanded = expand_n(2 * n - 1);
    u32 *a_expanded = expand_a(a, n, n_expanded);
    u32 *b_expanded = expand_a(b, n, n_expanded);

    // bit-reverse排列
    bit_reverse_permute(a_expanded, n_expanded);
    bit_reverse_permute(b_expanded, n_expanded);

    // 转换到计算域
    u32 omega_domain = mod_impl->to_compute_domain(omega);
    u32 w_domain = mod_impl->to_compute_domain(1);
    mod_impl->array_to_compute_domain(a_expanded, n_expanded);
    mod_impl->array_to_compute_domain(b_expanded, n_expanded);

    // 前向NTT
    ntt_forward_unified_omp(a_expanded, n_expanded, p, omega_domain, w_domain, mod_impl);
    ntt_forward_unified_omp(b_expanded, n_expanded, p, omega_domain, w_domain, mod_impl);

// 点乘（数据已经在计算域中）
#pragma omp parallel for
    for (u32 i = 0; i < n_expanded; ++i)
        a_expanded[i] = mod_impl->mul(a_expanded[i], b_expanded[i]);

    // 逆向NTT
    ntt_inverse_unified_omp(a_expanded, n_expanded, p, omega_domain, w_domain, mod_impl);

    // 转换回普通域
    mod_impl->array_from_compute_domain(a_expanded, n_expanded);

    // bit-reverse排列
    bit_reverse_permute(a_expanded, n_expanded);

    // 复制结果
    for (u32 i = 0; i < 2 * n - 1; ++i)
        ab[i] = a_expanded[i];

    // 清理内存
    delete[] a_expanded;
    delete[] b_expanded;
}
