#pragma once

#include "op.h"
#include "transform.h"
#include "type.h"

#define OMEGA 3 // 998244353 的原根

/**
 * @brief 使用NTT优化的多项式乘法
 *
 * 进一步实现了SIMD优化。
 *
 * @param a 多项式系数
 * @param b 多项式系数
 * @param ab 结果
 * @param n 多项式长度
 * @param p 模数（质数）
 */
inline void poly_multiply_ntt_simd(int *a, int *b, int *ab, int n, int p)
{
  Mod mod(p);

  u32 n_expanded = expand_n(2 * n - 1);
  u32 *a_expanded = expand_a((u32 *)a, n, n_expanded);
  u32 *b_expanded = expand_a((u32 *)b, n, n_expanded);

  ntt_forward(a_expanded, n_expanded, p, OMEGA);
  ntt_forward(b_expanded, n_expanded, p, OMEGA);
  for (u32 i = 0; i < n_expanded; ++i)
    ab[i] = mod.mul(a_expanded[i], b_expanded[i]);
  ntt_inverse((u32 *)ab, n_expanded, p, OMEGA);
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
inline void poly_multiply_ntt(int *a, int *b, int *ab, int n, int p)
{
  Mod mod(p);

  u32 n_expanded = expand_n(2 * n - 1);
  u32 *a_expanded = expand_a((u32 *)a, n, n_expanded);
  u32 *b_expanded = expand_a((u32 *)b, n, n_expanded);

  ntt_forward(a_expanded, n_expanded, p, OMEGA);
  ntt_forward(b_expanded, n_expanded, p, OMEGA);
  for (u32 i = 0; i < n_expanded; ++i)
    ab[i] = mod.mul(a_expanded[i], b_expanded[i]);
  ntt_inverse((u32 *)ab, n_expanded, p, OMEGA);
}

inline void poly_multiply_naive(int *a, int *b, int *ab, int n, int p)
{
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      ab[i + j] = (1LL * a[i] * b[j] % p + ab[i + j]) % p;
}