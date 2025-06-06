#pragma once

#include "../general/utils.h"
#include "op.h"

/**
 * @brief NTT 正变换：a(x) → A(ω)
 *
 * 输入的顺序是 bit-reversed，输出的顺序是自然顺序。
 *
 * 进行就地变换是缓存友好的。
 *
 * @param a 多项式系数，变换后表示频域系数
 * @param n 多项式长度（普通整数）
 * @param p 模数（普通整数）
 * @param omega 原根
 */
template <typename T>
inline void ntt_forward_Barrett(T *a, T n, T p, T omega)
{
    BarrettMod<T> mod(p, 100);

    for (T mid = 1; mid < n; mid <<= 1)
    {
        T Wn = mod.pow(omega, (p - 1) / (mid << 1));
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
    }
}

/**
 * @brief NTT 逆变换：A(ω) → a(x)
 *
 * 输入的顺序是自然顺序，输出的顺序是 bit-reversed。
 *
 * @param a 频域系数，变换后表示多项式系数
 * @param n 多项式长度
 * @param p 模数
 * @param omega 原根，已经是正变换的 ω，在调用时传 mod.inv(omega)
 */
template <typename T>
inline void ntt_inverse_Barrett(T *a, T n, T p, T omega)
{
    BarrettMod<T> mod(p, 100);

    for (T mid = n >> 1; mid > 0; mid >>= 1)
    {
        T Wn = mod.pow(omega, (p - 1) / (mid << 1)); // Wn = ω⁻¹^((p-1)/(2*mid))
        for (T j = 0; j < n; j += (mid << 1))
        {
            T w = 1;
            for (T k = 0; k < mid; ++k, w = mod.mul(w, Wn))
            {
                T x = a[j + k];
                T y = a[j + k + mid];
                a[j + k] = mod.add(x, y);
                a[j + k + mid] = mod.mul(w, mod.sub(x, y));
            }
        }
    }

    T inv_n = mod.inv(n);
    for (T i = 0; i < n; ++i)
        a[i] = mod.mul(a[i], inv_n);
}

/**
 * @brief 使用NTT优化的多项式乘法
 *
 * @param a 多项式系数
 * @param b 多项式系数
 * @param ab 结果
 * @param n 多项式长度
 * @param p 模数（质数）
 */
template <typename T>
inline void poly_multiply_ntt_Barrett(T *a, T *b, T *ab, T n, T p, T OMEGA = 3)
{
    BarrettMod<T> mod(p, 100);

    T n_expanded = expand_n(2 * n - 1);
    T *a_expanded = expand_a((T *)a, n, n_expanded);
    T *b_expanded = expand_a((T *)b, n, n_expanded);

    bit_reverse_permute(a_expanded, n_expanded);
    bit_reverse_permute(b_expanded, n_expanded);

    T *a_copy = new T[n_expanded]{};
    T *b_copy = new T[n_expanded]{};
    T *ab_result = new T[n_expanded]{};

    for (T i = 0; i < n_expanded; ++i)
        a_copy[i] = a_expanded[i];
    for (T i = 0; i < n_expanded; ++i)
        b_copy[i] = b_expanded[i];

    ntt_forward_Barrett(a_copy, n_expanded, p, OMEGA);
    ntt_forward_Barrett(b_copy, n_expanded, p, OMEGA);

    for (T i = 0; i < n_expanded; ++i)
        ab_result[i] = mod.mul(a_copy[i], b_copy[i]);

    ntt_inverse_Barrett(ab_result, n_expanded, p, mod.inv(OMEGA));

    for (T i = 0; i < n_expanded; ++i)
        ab[i] = ab_result[i];

    bit_reverse_permute((T *)ab, n_expanded);

    delete[] a_copy;
    delete[] b_copy;
    delete[] ab_result;
}