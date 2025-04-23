#pragma once

#include "op.h"
#include "transform.h"

#define PRIMITIVE_ROOT 3 // 998244353 的原根

/**
 * @brief 输入两个多项式系数（模 mod）
 * 返回多项式卷积结果（模 mod）
 * 自动补零至 2 的幂，执行 NTT、点乘、逆NTT
 * 内部调用顺序：
 * - root 表构造
 * - NTT 变换（正）
 * - 点乘（可并行）
 * - NTT 变换（逆）
 * - 缩减结果长度
 * @param a 多项式系数
 * @param b 多项式系数
 * @param ab 结果
 * @param n 多项式长度
 * @param p 模数（质数）
 */
void poly_multiply_ntt(int *a, int *b, int *ab, int n, int p) {
  ntt_forward(a, n, p, PRIMITIVE_ROOT);
  ntt_forward(b, n, p, PRIMITIVE_ROOT);
  for (int i = 0; i < n; ++i) {
    ab[i] = mod_mul(a[i], b[i], p);
  }
  ntt_inverse(ab, n, p, PRIMITIVE_ROOT);
}

void poly_multiply_naive(int *a, int *b, int *ab, int n, int p) {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      ab[i + j] = (1LL * a[i] * b[j] % p + ab[i + j]) % p;
    }
  }
}