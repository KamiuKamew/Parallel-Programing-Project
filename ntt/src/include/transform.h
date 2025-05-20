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
template <typename T>
inline void ntt_forward(T *a, T n, T p, T omega)
{
  Mod<T> mod(p);

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
 * @param a 频域系数，变换后表示多项式系数
 * @param n 多项式长度
 * @param p 模数
 * @param omega 原根
 */
template <typename T>
inline void ntt_inverse(T *a, T n, T p, T omega)
{
  Mod<T> mod(p);

  ntt_forward(a, n, p, mod.inv(omega));

  for (T i = 0; i < n; ++i)
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
template <typename T>
inline void ntt_forward_mont(T *a_mont, T n, T p, T omega_mont)
{
  using T_mont = T;

  MontMod<T> montMod(p);

  for (T mid = 1; mid < n; mid <<= 1)
  {
    T_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1));
    for (T j = 0; j < n; j += (mid << 1))
    {
      T_mont w_mont = montMod.from_T(1);
      for (T k = 0; k < mid; ++k, w_mont = montMod.mul(w_mont, Wn_mont))
      {
        T_mont x_mont = a_mont[j + k];
        T_mont y_mont = montMod.mul(w_mont, a_mont[j + k + mid]);
        a_mont[j + k] = montMod.add(x_mont, y_mont);
        a_mont[j + k + mid] = montMod.sub(x_mont, y_mont);
      }
    }
  }
}

template <typename T>
inline void ntt_forward_mont_before_simd(T *a_mont, T n, T p, T omega_mont)
{
  using T_mont = T;

  MontMod<T> montMod(p);

  for (T mid = 1; mid < n; mid <<= 1)
  {
    switch (mid)
    {
    case 1:
    {
      for (T j = 0; j < n; j += (mid << 1))
      {
        T_mont w_mont = montMod.from_T(1);
        T_mont x_mont = a_mont[j];
        T_mont y_mont = montMod.mul(w_mont, a_mont[j + 1]);
        a_mont[j] = montMod.add(x_mont, y_mont);
        a_mont[j + 1] = montMod.sub(x_mont, y_mont);
      }
      break;
    }

    case 2:
    {
      T_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1));
      for (T j = 0; j < n; j += (mid << 1))
      {
        T_mont w_mont_0 = montMod.from_T(1);
        T_mont w_mont_1 = montMod.mul(w_mont_0, Wn_mont);
        T_mont x_mont_0 = a_mont[j + 0];
        T_mont x_mont_1 = a_mont[j + 1];
        T_mont y_mont_0 = montMod.mul(w_mont_0, a_mont[j + 2]);
        T_mont y_mont_1 = montMod.mul(w_mont_1, a_mont[j + 3]);
        a_mont[j + 0] = montMod.add(x_mont_0, y_mont_0);
        a_mont[j + 1] = montMod.add(x_mont_1, y_mont_1);
        a_mont[j + 2] = montMod.sub(x_mont_0, y_mont_0);
        a_mont[j + 3] = montMod.sub(x_mont_1, y_mont_1);
      }
      break;
    }

    default: // mid >= 4, parallelizable
    {
      T_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1));
      for (T j = 0; j < n; j += (mid << 1))
      {
        T_mont w_mont_0 = montMod.from_T(1);
        T_mont w_mont_1 = montMod.mul(w_mont_0, Wn_mont);
        T_mont w_mont_2 = montMod.mul(w_mont_1, Wn_mont);
        T_mont w_mont_3 = montMod.mul(w_mont_2, Wn_mont);
        T_mont Wn_mont_4 = montMod.pow(Wn_mont, 4);
        for (T k = 0; k < mid; k += 4)
        {
          T_mont x_mont_0 = a_mont[j + k + 0];
          T_mont x_mont_1 = a_mont[j + k + 1];
          T_mont x_mont_2 = a_mont[j + k + 2];
          T_mont x_mont_3 = a_mont[j + k + 3];

          T_mont y_mont_0 = montMod.mul(w_mont_0, a_mont[j + k + mid + 0]);
          T_mont y_mont_1 = montMod.mul(w_mont_1, a_mont[j + k + mid + 1]);
          T_mont y_mont_2 = montMod.mul(w_mont_2, a_mont[j + k + mid + 2]);
          T_mont y_mont_3 = montMod.mul(w_mont_3, a_mont[j + k + mid + 3]);

          a_mont[j + k + 0] = montMod.add(x_mont_0, y_mont_0);
          a_mont[j + k + 1] = montMod.add(x_mont_1, y_mont_1);
          a_mont[j + k + 2] = montMod.add(x_mont_2, y_mont_2);
          a_mont[j + k + 3] = montMod.add(x_mont_3, y_mont_3);

          a_mont[j + k + mid + 0] = montMod.sub(x_mont_0, y_mont_0);
          a_mont[j + k + mid + 1] = montMod.sub(x_mont_1, y_mont_1);
          a_mont[j + k + mid + 2] = montMod.sub(x_mont_2, y_mont_2);
          a_mont[j + k + mid + 3] = montMod.sub(x_mont_3, y_mont_3);

          w_mont_0 = montMod.mul(w_mont_0, Wn_mont_4);
          w_mont_1 = montMod.mul(w_mont_1, Wn_mont_4);
          w_mont_2 = montMod.mul(w_mont_2, Wn_mont_4);
          w_mont_3 = montMod.mul(w_mont_3, Wn_mont_4);
        }
      }
      break;
    }
    }
  }
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
template <typename T>
inline void ntt_inverse_mont(T *a_mont, T n, T p, T omega_mont)
{
  using T_mont = T;

  MontMod<T> montMod(p);

  for (T mid = n >> 1; mid > 0; mid >>= 1)
  {
    T_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1)); // Wn = ω⁻¹^((p-1)/(2*mid))
    for (T j = 0; j < n; j += (mid << 1))
    {
      T_mont w_mont = montMod.from_T(1);
      for (T k = 0; k < mid; ++k, w_mont = montMod.mul(w_mont, Wn_mont))
      {
        T_mont x_mont = a_mont[j + k];
        T_mont y_mont = a_mont[j + k + mid];
        a_mont[j + k] = montMod.add(x_mont, y_mont);
        a_mont[j + k + mid] = montMod.mul(w_mont, montMod.sub(x_mont, y_mont));
      }
    }
  }

  T_mont inv_n = montMod.inv(montMod.from_T(n));
  for (T i = 0; i < n; ++i)
    a_mont[i] = montMod.mul(a_mont[i], inv_n);
}

template <typename T>
inline void ntt_inverse_mont_before_simd(T *a_mont, T n, T p, T omega_mont)
{
  using T_mont = T;

  MontMod<T> montMod(p);

  for (T mid = n >> 1; mid > 0; mid >>= 1)
  {
    T_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1)); // Wn = ω⁻¹^((p-1)/(2*mid))
    switch (mid)
    {
    case 1:
    {
      for (T j = 0; j < n; j += (mid << 1))
      {
        T_mont w_mont = montMod.from_T(1);
        T_mont x_mont = a_mont[j + 0];
        T_mont y_mont = a_mont[j + 1];
        a_mont[j + 0] = montMod.add(x_mont, y_mont);
        a_mont[j + 1] = montMod.mul(w_mont, montMod.sub(x_mont, y_mont));
      }
      break;
    }

    case 2:
    {
      for (T j = 0; j < n; j += (mid << 1))
      {
        T_mont w_mont_0 = montMod.from_T(1);
        T_mont w_mont_1 = montMod.mul(w_mont_0, Wn_mont);
        T_mont x_mont_0 = a_mont[j + 0];
        T_mont x_mont_1 = a_mont[j + 1];
        T_mont y_mont_0 = a_mont[j + 2];
        T_mont y_mont_1 = a_mont[j + 3];
        a_mont[j + 0] = montMod.add(x_mont_0, y_mont_0);
        a_mont[j + 1] = montMod.add(x_mont_1, y_mont_1);
        a_mont[j + 2] = montMod.mul(w_mont_0, montMod.sub(x_mont_0, y_mont_0));
        a_mont[j + 3] = montMod.mul(w_mont_1, montMod.sub(x_mont_1, y_mont_1));
      }
      break;
    }

    default: // mid >= 4, parallelizable
    {
      for (T j = 0; j < n; j += (mid << 1))
      {
        T_mont w_mont_0 = montMod.from_T(1);
        T_mont w_mont_1 = montMod.mul(w_mont_0, Wn_mont);
        T_mont w_mont_2 = montMod.mul(w_mont_1, Wn_mont);
        T_mont w_mont_3 = montMod.mul(w_mont_2, Wn_mont);
        T_mont Wn_mont_4 = montMod.pow(Wn_mont, 4);
        for (T k = 0; k < mid; k += 4)
        {
          T_mont x_mont_0 = a_mont[j + k + 0];
          T_mont x_mont_1 = a_mont[j + k + 1];
          T_mont x_mont_2 = a_mont[j + k + 2];
          T_mont x_mont_3 = a_mont[j + k + 3];

          T_mont y_mont_0 = a_mont[j + k + mid + 0];
          T_mont y_mont_1 = a_mont[j + k + mid + 1];
          T_mont y_mont_2 = a_mont[j + k + mid + 2];
          T_mont y_mont_3 = a_mont[j + k + mid + 3];

          a_mont[j + k + 0] = montMod.add(x_mont_0, y_mont_0);
          a_mont[j + k + 1] = montMod.add(x_mont_1, y_mont_1);
          a_mont[j + k + 2] = montMod.add(x_mont_2, y_mont_2);
          a_mont[j + k + 3] = montMod.add(x_mont_3, y_mont_3);

          a_mont[j + k + mid + 0] = montMod.mul(w_mont_0, montMod.sub(x_mont_0, y_mont_0));
          a_mont[j + k + mid + 1] = montMod.mul(w_mont_1, montMod.sub(x_mont_1, y_mont_1));
          a_mont[j + k + mid + 2] = montMod.mul(w_mont_2, montMod.sub(x_mont_2, y_mont_2));
          a_mont[j + k + mid + 3] = montMod.mul(w_mont_3, montMod.sub(x_mont_3, y_mont_3));

          w_mont_0 = montMod.mul(w_mont_0, Wn_mont_4);
          w_mont_1 = montMod.mul(w_mont_1, Wn_mont_4);
          w_mont_2 = montMod.mul(w_mont_2, Wn_mont_4);
          w_mont_3 = montMod.mul(w_mont_3, Wn_mont_4);
        }
      }
      break;
    }
    }
  }

  T_mont inv_n = montMod.inv(montMod.from_T(n));
  for (T i = 0; i < n; ++i)
    a_mont[i] = montMod.mul(a_mont[i], inv_n);
}
