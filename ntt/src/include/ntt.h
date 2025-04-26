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
  MontMod montMod(p);

  u32 n_expanded = expand_n(2 * n - 1);
  u32 *a_expanded = expand_a((u32 *)a, n, n_expanded);
  u32 *b_expanded = expand_a((u32 *)b, n, n_expanded);

  u32_mont *a_mont = new u32_mont[n_expanded];
  u32_mont *b_mont = new u32_mont[n_expanded];
  u32_mont *ab_mont = new u32_mont[n_expanded];
  for (u32 i = 0; i < n_expanded; ++i)
    a_mont[i] = montMod.from_u32(a_expanded[i]);
  for (u32 i = 0; i < n_expanded; ++i)
    b_mont[i] = montMod.from_u32(b_expanded[i]);
  u32_mont omega_mont = montMod.from_u32(OMEGA);

  bit_reverse_permute(a_mont, n_expanded);
  bit_reverse_permute(b_mont, n_expanded);
  ntt_forward_mont(a_mont, n_expanded, p, omega_mont);
  ntt_forward_mont(b_mont, n_expanded, p, omega_mont);

  for (u32 i = 0; i < n_expanded; ++i)
    ab_mont[i] = montMod.mul(a_mont[i], b_mont[i]);

  ntt_inverse_dit_mont(ab_mont, n_expanded, p, montMod.inv(omega_mont));
  bit_reverse_permute(ab_mont, n_expanded);

  for (u32 i = 0; i < n_expanded; ++i)
    ab[i] = montMod.to_u32(ab_mont[i]);
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
  MontMod montMod(p);

  u32 n_expanded = expand_n(2 * n - 1);
  u32 *a_expanded = expand_a((u32 *)a, n, n_expanded);
  u32 *b_expanded = expand_a((u32 *)b, n, n_expanded);

  u32_mont *a_mont = new u32_mont[n_expanded];
  u32_mont *b_mont = new u32_mont[n_expanded];
  u32_mont *ab_mont = new u32_mont[n_expanded];
  for (u32 i = 0; i < n_expanded; ++i)
    a_mont[i] = montMod.from_u32(a_expanded[i]);
  for (u32 i = 0; i < n_expanded; ++i)
    b_mont[i] = montMod.from_u32(b_expanded[i]);
  u32_mont omega_mont = montMod.from_u32(OMEGA);

  bit_reverse_permute(a_mont, n_expanded);
  bit_reverse_permute(b_mont, n_expanded);
  ntt_forward_mont(a_mont, n_expanded, p, omega_mont);
  ntt_forward_mont(b_mont, n_expanded, p, omega_mont);

  for (u32 i = 0; i < n_expanded; ++i)
    ab_mont[i] = montMod.mul(a_mont[i], b_mont[i]);

  ntt_inverse_dit_mont(ab_mont, n_expanded, p, montMod.inv(omega_mont));
  bit_reverse_permute(ab_mont, n_expanded);

  for (u32 i = 0; i < n_expanded; ++i)
    ab[i] = montMod.to_u32(ab_mont[i]);
}

inline void poly_multiply_naive(int *a, int *b, int *ab, int n, int p)
{
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      ab[i + j] = (1LL * a[i] * b[j] % p + ab[i + j]) % p;
}