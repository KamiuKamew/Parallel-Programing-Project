#pragma once

#include "../general/op.h"
#include "../general/utils.h"
#include <algorithm>
#include <iostream>

/** 将 n 扩展为 4 的幂。 */
template <typename T>
inline T expand_n_radix4(T n)
{
  T res = 1;
  while (res < n)
    res <<= 2; // Multiply by 4
  return res;
}

/** 将 a 的长度扩展为 4 的幂。 */
template <typename T>
inline T *expand_a_radix4(T *a, T current_len, T expanded_len)
{
  T *a_expanded = new T[expanded_len];
  for (T i = 0; i < current_len; ++i)
    a_expanded[i] = a[i];
  for (T i = current_len; i < expanded_len; ++i)
    a_expanded[i] = 0;
  return a_expanded;
}

/**
 * @brief 基-4 NTT 正变换：a(x) → A(ω)
 *
 * @param a 多项式系数，变换后表示频域系数
 * @param n 多项式长度
 * @param p 模数
 * @param omega 原根
 */
template <typename T>
inline void ntt_forward_radix4(T *a, T n, T p, T omega)
{
  Mod<T> mod(p);

  for (T len = 4; len <= n; len <<= 2)
  {
    T Wn = mod.pow(omega, (p - 1) / len);
    for (T i = 0; i < n; i += len)
    {
      T w = 1;
      for (T j = 0; j < len / 4; ++j)
      {
        T t0 = a[i + j];
        T t1 = mod.mul(a[i + j + len / 4], w);
        T t2 = mod.mul(a[i + j + len / 2], mod.mul(w, w));
        T t3 = mod.mul(a[i + j + 3 * len / 4], mod.mul(mod.mul(w, w), w));

        T A = mod.add(t0, t2);
        T B = mod.sub(t0, t2);
        T C = mod.add(t1, t3);
        T D = mod.sub(t1, t3);

        a[i + j] = mod.add(A, C);
        a[i + j + len / 4] = mod.add(B, D);
        a[i + j + len / 2] = mod.sub(A, C);
        a[i + j + 3 * len / 4] = mod.sub(B, D);

        w = mod.mul(w, Wn);
      }
    }
  }
}

/**
 * @brief 基-4 NTT 逆变换：A(ω) → a(x)
 *
 * @param a 频域系数，变换后表示多项式系数
 * @param n 多项式长度
 * @param p 模数
 * @param omega 原根
 */
template <typename T>
void ntt_inverse_radix4(T *a, T n, T p, T omega)
{
  Mod<T> mod(p);
  T omega_inv = mod.inv(omega);

  for (T len = n; len >= 4; len >>= 2)
  {
    T Wn_inv = mod.pow(omega_inv, (p - 1) / len);
    for (T i = 0; i < n; i += len)
    {
      T w = 1;
      for (T j = 0; j < len / 4; ++j)
      {
        T t0 = a[i + j];
        T t1 = a[i + j + len / 4];
        T t2 = a[i + j + len / 2];
        T t3 = a[i + j + 3 * len / 4];

        T A = mod.add(t0, t2);
        T B = mod.sub(t0, t2);
        T C = mod.add(t1, t3);
        T D = mod.sub(t1, t3);

        a[i + j] = mod.add(A, C);
        a[i + j + len / 4] = mod.mul(mod.add(B, D), w);
        a[i + j + len / 2] = mod.mul(mod.sub(A, C), mod.mul(w, w));
        a[i + j + 3 * len / 4] = mod.mul(mod.sub(B, D), mod.mul(mod.mul(w, w), w));

        w = mod.mul(w, Wn_inv);
      }
    }
  }

  T inv_n = mod.inv(n);
  for (T i = 0; i < n; ++i)
    a[i] = mod.mul(a[i], inv_n);
}

/**
 * @brief 使用基-4 NTT优化的多项式乘法
 *
 * @param a 多项式系数
 * @param b 多项式系数
 * @param ab 结果
 * @param n_a 多项式a的长度
 * @param n_b 多项式b的长度
 * @param p 模数（质数）
 */
template <typename T>
inline void poly_multiply_ntt_radix4(T *a, T *b, T *ab, T n_a, T n_b, T p, T OMEGA = 3)
{
  Mod<T> mod(p);

  T n_result_len = n_a + n_b - 1;               // The actual length of the result polynomial
  T n_expanded = expand_n_radix4(n_result_len); // Smallest power of 4 >= n_result_len

  T *a_expanded = expand_a_radix4((T *)a, n_a, n_expanded);
  T *b_expanded = expand_a_radix4((T *)b, n_b, n_expanded);

  // DIT NTT requires bit-reversed input
  bit_reverse_permute(a_expanded, n_expanded);
  bit_reverse_permute(b_expanded, n_expanded);

  T *ab_expanded = new T[n_expanded]{};

  ntt_forward_radix4(a_expanded, n_expanded, p, OMEGA);
  ntt_forward_radix4(b_expanded, n_expanded, p, OMEGA);

  for (T i = 0; i < n_expanded; ++i)
    ab_expanded[i] = mod.mul(a_expanded[i], b_expanded[i]);

  ntt_inverse_radix4(ab_expanded, n_expanded, p, OMEGA);

  // Apply bit reversal to get correct order
  bit_reverse_permute(ab_expanded, n_expanded);

  // Copy results back to ab, only up to n_result_len
  for (T i = 0; i < n_result_len; ++i)
    ab[i] = ab_expanded[i];

  delete[] a_expanded;
  delete[] b_expanded;
  delete[] ab_expanded;
}
