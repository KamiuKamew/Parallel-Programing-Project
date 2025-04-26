#pragma once

#include "general/op.h"

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
 * 输入的顺序是 bit-reversed，输出的顺序是自然顺序。
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
 * 复用正变换的版本。
 * 输入的顺序是 bit-reversed，输出的顺序是自然顺序。
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

/**
 * @brief NTT 逆变换：A(ω) → a_mont(x)
 *
 * 输入的顺序是自然顺序，输出的顺序是 bit-reversed。
 *
 * @param a_mont 频域系数，变换后表示多项式系数
 * @param n 多项式长度
 * @param p 模数
 * @param omega_mont 原根，已经是正变换的 ω，在调用时传 montMod.inv(omega_mont)
 */
inline void ntt_inverse_dit_mont(u32_mont *a_mont, u32 n, u32 p, u32_mont omega_mont)
{
  MontMod montMod(p);

  for (u32 mid = n >> 1; mid > 0; mid >>= 1)
  {
    u32_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1)); // Wn = ω⁻¹^((p-1)/(2*mid))
    for (u32 j = 0; j < n; j += (mid << 1))
    {
      u32_mont w_mont = montMod.from_u32(1);
      for (u32 k = 0; k < mid; ++k, w_mont = montMod.mul(w_mont, Wn_mont))
      {
        u32_mont x_mont = a_mont[j + k];
        u32_mont y_mont = a_mont[j + k + mid];
        a_mont[j + k] = montMod.add(x_mont, y_mont);
        a_mont[j + k + mid] = montMod.mul(w_mont, montMod.sub(x_mont, y_mont));
      }
    }
  }

  u32_mont inv_n = montMod.inv(montMod.from_u32(n));
  for (u32 i = 0; i < n; ++i)
    a_mont[i] = montMod.mul(a_mont[i], inv_n);
}
