#pragma once

#include "type.h"
#include "op.h"
#include "transform.h"

#define OMEGA 3 // 998244353 的原根

/**
 * @brief 使用NTT优化的多项式乘法
 * 自动补零至 2 的幂，执行 NTT、点乘、逆NTT
 * @param a 多项式系数
 * @param b 多项式系数
 * @param ab 结果
 * @param n 多项式长度
 * @param p 模数（质数）
 */
void poly_multiply_ntt(int *a, int *b, int *ab, int n, int p) {
  u32 n_expanded = expand_n(2 * n - 1);
  u32 *a_expanded = expand_a((u32*)a, n, n_expanded);
  u32 *b_expanded = expand_a((u32*)b, n, n_expanded);

  ntt_forward(a_expanded, n_expanded, p, OMEGA);
  ntt_forward(b_expanded, n_expanded, p, OMEGA);
  for (int i = 0; i < n_expanded; ++i) {
    ab[i] = mod_mul(a_expanded[i], b_expanded[i], p);
  }
  ntt_inverse((u32*)ab, n_expanded, p, OMEGA);
}

void poly_multiply_naive(int *a, int *b, int *ab, int n, int p) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      ab[i + j] = (1LL * a[i] * b[j] % p + ab[i + j]) % p;
    }
  }
}