#pragma once

#include "op.h"

/**
 * @brief 就地对 a 做 bit-reverse 置换
 *
 * @param a 输入序列
 * @param n 序列长度
 */
void bit_reverse_permute(int *a, int n) {
  int lg_n = 0;
  while ((1 << lg_n) < n)
    ++lg_n;

  for (int i = 0; i < n; ++i) {
    int j = 0;
    for (int k = 0; k < lg_n; ++k) {
      if (i & (1 << k)) {
        j |= (1 << (lg_n - 1 - k));
      }
    }
    if (i < j) {
      auto tmp = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
  }
}

/**
 * @brief NTT 正变换：a(x) → A(ω)
 *
 * 进行就地变换是缓存友好的。
 *
 * @param a 多项式系数，变换后表示频域系数
 * @param n 多项式长度
 * @param p 模数
 * @param omega 原根
 */
void ntt_forward(int *a, int n, int p, int omega) {
  bit_reverse_permute(a, n); // 注意 n 一定是 2 的幂

  for (int mid = 1; mid < n; mid <<= 1) {
    int Wn = mod_pow(omega, (p - 1) / (mid << 1), p);
    for (int j = 0; j < n; j += (mid << 1)) {
      int w = 1;
      for (int k = 0; k < mid; ++k, w = mod_mul(w, Wn, p)) {
        int x = a[j + k];
        int y = mod_mul(w, a[j + k + mid], p);
        a[j + k] = mod_add(x, y, p);
        a[j + k + mid] = mod_sub(x, y, p);
      }
    }
  }
}

/**
 * @brief NTT 逆变换：A(ω) → a(x)
 *
 * @param a 频域系数，变换后表示多项式系数
 * @param n 多项式长度
 * @param p 模数
 * @param omega 原根
 */
void ntt_inverse(int *a, int n, int p, int omega) {
  bit_reverse_permute(a, n);

  int inv_root = mod_inv(omega, p);
  for (int mid = 1; mid < n; mid <<= 1) {
    int Wn = mod_pow(inv_root, (p - 1) / (mid << 1), p);
    for (int j = 0; j < n; j += (mid << 1)) {
      int w = 1;
      for (int k = 0; k < mid; ++k, w = mod_mul(w, Wn, p)) {
        int x = a[j + k];
        int y = mod_mul(w, a[j + k + mid], p);
        a[j + k] = mod_add(x, y, p);
        a[j + k + mid] = mod_sub(x, y, p);
      }
    }
  }

  // 最后每个元素乘以 n 的逆元
  int inv_n = mod_inv(n, p);
  for (int i = 0; i < n; ++i) {
    a[i] = mod_mul(a[i], inv_n, p);
  }
}
