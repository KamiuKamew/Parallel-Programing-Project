#pragma once

#include "op.h"
#include "type.h"

/** 将 n 扩展为 2 的幂。 */
inline u32 expand_n(u32 n)
{
  u32 lg_n = 0;
  while ((1u << lg_n) < n)
    ++lg_n;
  return 1 << lg_n;
}

/** 将 a 的长度扩展为 2 的幂。 */
inline u32 *expand_a(u32 *a, u32 n, u32 n_expanded)
{
  // 虽然但是，这里应该有一个检查n是不是2的幂的逻辑
  // 话又说回来，这样得用别的库
  // 然而我又不想用别的库
  // 所以这里就先不检查了

  u32 *a_expanded = new u32[n_expanded];
  for (u32 i = 0; i < n; ++i)
    a_expanded[i] = a[i];
  for (u32 i = n; i < n_expanded; ++i)
    a_expanded[i] = 0;
  return a_expanded;
}

/** 将 a 转换成 SIMD 类型。 */
inline u32x4 *to_simd(u32 *a, u32 n_expanded)
{
  u32x4 *a_simd = new u32x4[n_expanded / 4];
  for (u32 i = 0; i < n_expanded / 4; ++i)
    a_simd[i] = vld1q_u32(&a[i * 4]);
  return a_simd;
}

/** 将 a_simd 转换成普通类型。 */
inline u32 *from_simd(u32x4 *a_simd, u32 n_expanded)
{
  u32 *a = new u32[n_expanded];
  for (u32 i = 0; i < n_expanded / 4; ++i)
    vst1q_u32(&a[i * 4], a_simd[i]);
  return a;
}

/**
 * @brief 就地对 a 做 bit-reverse 置换。
 *
 * 这个函数看起来不像是能 SIMD 优化的。
 * 而且这个函数应该在向量化之前用。
 *
 * @param a 输入序列（是不是 Montgomery 数域无所谓）
 * @param n 序列长度（扩展过，是2的幂）
 */
inline void bit_reverse_permute(u32 *a, u32 n)
{
  u32 lg_n = 0;
  while ((1u << lg_n) < n)
    ++lg_n;

  for (u32 i = 0; i < n; ++i)
  {
    u32 j = 0;
    for (u32 k = 0; k < lg_n; ++k)
      if (i & (1 << k))
        j |= (1 << (lg_n - 1 - k));
    if (i < j)
    {
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
inline void ntt_forward(u32 *a, u32 n, u32 p, u32 omega)
{
  Mod mod(p);

  for (u32 mid = 1; mid < n; mid <<= 1)
  {
    u32 Wn = mod.pow(omega, (p - 1) / (mid << 1));
    for (u32 j = 0; j < n; j += (mid << 1))
    {
      u32 w = 1;
      for (u32 k = 0; k < mid; ++k, w = mod.mul(w, Wn))
      {
        u32 x = a[j + k];
        u32 y = mod.mul(w, a[j + k + mid]);
        a[j + k] = mod.add(x, y);
        a[j + k + mid] = mod.sub(x, y);
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
inline void ntt_inverse(u32 *a, u32 n, u32 p, u32 omega)
{
  Mod mod(p);

  ntt_forward(a, n, p, mod.inv(omega));

  for (u32 i = 0; i < n; ++i)
    a[i] = mod.mul(a[i], mod.inv(n)); // 最后每个元素乘以 n 的逆元
}

/**
 * @brief NTT 正变换：a(x) → A(ω)
 *
 * 进行就地变换是缓存友好的。
 *
 * @param a_mont 多项式系数（位于 Montgomery 数域），变换后表示频域系数（位于 Montgomery 数域）
 * @param n 多项式长度（普通整数）
 * @param p 模数（普通整数）
 * @param omega_mont 原根（位于 Montgomery 数域）
 */
inline void ntt_forward_mont(u32_mont *a_mont, u32 n, u32 p, u32_mont omega_mont)
{
  MontMod montMod(p);

  for (u32 mid = 1; mid < n; mid <<= 1)
  {
    u32_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1));
    for (u32 j = 0; j < n; j += (mid << 1))
    {
      u32_mont w_mont = montMod.from_u32(1);
      for (u32 k = 0; k < mid; ++k, w_mont = montMod.mul(w_mont, Wn_mont))
      {
        u32_mont x_mont = a_mont[j + k];
        u32_mont y_mont = montMod.mul(w_mont, a_mont[j + k + mid]);
        a_mont[j + k] = montMod.add(x_mont, y_mont);
        a_mont[j + k + mid] = montMod.sub(x_mont, y_mont);
      }
    }
  }
}

/**
 * @brief NTT 逆变换：A(ω) → a_mont(x)
 *
 * @param a_mont 频域系数，变换后表示多项式系数
 * @param n 多项式长度
 * @param p 模数
 * @param omega_mont 原根
 */
inline void ntt_inverse_mont(u32_mont *a_mont, u32 n, u32 p, u32_mont omega_mont)
{
  MontMod montMod(p);

  ntt_forward_mont(a_mont, n, p, montMod.inv(omega_mont));

  u32_mont n_mont = montMod.from_u32(n);
  u32_mont inv_n_mont = montMod.inv(n_mont);
  for (u32 i = 0; i < n; ++i)
    a_mont[i] = montMod.mul(a_mont[i], inv_n_mont); // 最后每个元素乘以 n 的逆元
}