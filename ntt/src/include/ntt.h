#pragma once

#include "general/utils.h"
#include "transform.h"

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
inline void poly_multiply_ntt(T *a, T *b, T *ab, T n, T p, T OMEGA = 3)
{
  using T_mont = T;

  MontMod<T> montMod(p);

  T n_expanded = expand_n(2 * n - 1);
  T *a_expanded = expand_a((T *)a, n, n_expanded);
  T *b_expanded = expand_a((T *)b, n, n_expanded);

  bit_reverse_permute(a_expanded, n_expanded);
  bit_reverse_permute(b_expanded, n_expanded);

  T_mont *a_mont = new T_mont[n_expanded]{};
  T_mont *b_mont = new T_mont[n_expanded]{};
  T_mont *ab_mont = new T_mont[n_expanded]{};
  for (T i = 0; i < n_expanded; ++i)
    a_mont[i] = montMod.from_T(a_expanded[i]);
  for (T i = 0; i < n_expanded; ++i)
    b_mont[i] = montMod.from_T(b_expanded[i]);
  T_mont omega_mont = montMod.from_T(OMEGA);

  ntt_forward_mont(a_mont, n_expanded, p, omega_mont);
  ntt_forward_mont(b_mont, n_expanded, p, omega_mont);

  for (T i = 0; i < n_expanded; ++i)
    ab_mont[i] = montMod.mul(a_mont[i], b_mont[i]);

  ntt_inverse_mont(ab_mont, n_expanded, p, montMod.inv(omega_mont));

  for (T i = 0; i < n_expanded; ++i)
    ab[i] = montMod.to_T(ab_mont[i]);

  bit_reverse_permute((T *)ab, n_expanded);
}

inline void poly_multiply_naive(int *a, int *b, int *ab, int n, int p)
{
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      ab[i + j] = (1LL * a[i] * b[j] % p + ab[i + j]) % p;
}