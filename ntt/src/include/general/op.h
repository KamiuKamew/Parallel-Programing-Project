#pragma once

#include "type.h"

// === 普通模运算 ===
template <typename T>
class Mod
{
  using T2 = t_widen<T>;

public:
  Mod(T _mod) : mod(_mod) {}
  Mod(const Mod &) = delete;
  Mod &operator=(const Mod &) = delete;

  T add(T a, T b) const { return (a + b) % mod; }
  T sub(T a, T b) const { return (a >= b) ? (a - b) : (a + mod - b); }
  T mul(T a, T b) const { return ((T2)a * (T2)b) % mod; } // TODO：当 T = u128 时，a * b 可能溢出
  T pow(T base, T exp) const
  {
    T result = 1;
    while (exp > 0)
    {
      if (exp & 1)
      {
        result = mul(result, base);
      }
      base = mul(base, base);
      exp >>= 1;
    }
    return result;
  }
  T inv(T x) const { return pow(x, mod - 2); }

private:
  T mod;
};

using Mod128 = Mod<u128>;

// === Montgomery 模运算 ===
template <typename T>
class MontMod
{
  using T_mont = T;
  using T2 = t_widen<T>;
  static constexpr int word_bits = sizeof(T) * 8;

public:
  MontMod(T _mod) : mod(_mod)
  {
    // 计算 R^2 mod mod，这里 R = 2^mont_word_bits
    // First, calculate R_mod_n = 2^mont_word_bits mod mod
    T2 r_val_mod_n = 1;
    for (int i = 0; i < word_bits; ++i)
      r_val_mod_n = (r_val_mod_n << 1) % mod;
    // Now r_val_mod_n = 2^mont_word_bits mod mod
    // Then r2 = (2^mont_word_bits mod mod)^2 mod mod = (2^mont_word_bits)^2 mod mod
    r2 = (T2)r_val_mod_n * r_val_mod_n % mod;

    // 计算 -mod^(-1) mod R (R = 2^mont_word_bits)
    // 即满足 R*R_inv - mod*m_inv = 1 (mod R)
    T inv = 1;
    // 牛顿迭代法求模逆 mod R
    int inv_iterations;
    if (word_bits <= 2)
      inv_iterations = 1;
    else if (word_bits <= 4)
      inv_iterations = 2;
    else if (word_bits <= 8)
      inv_iterations = 3;
    else if (word_bits <= 16)
      inv_iterations = 4;
    else if (word_bits <= 32)
      inv_iterations = 5;
    else if (word_bits <= 64)
      inv_iterations = 6;
    else if (word_bits <= 128)
      inv_iterations = 7;
    else
      inv_iterations = 5;

    for (int i = 0; i < inv_iterations; ++i)
      inv = (T)((T2)inv * (2 - (T2)mod * inv));
    neg_r_inv = -inv;
  }

  MontMod(const MontMod &) = delete;
  MontMod &operator=(const MontMod &) = delete;

  T_mont from_T(T a) const { return reduce((T2)a * r2); }
  T to_T(T_mont a_mont) const { return reduce((T2)a_mont); }

  // Montgomery规约，计算 t * R^(-1) mod mod, where R = 2^mont_word_bits
  T_mont reduce(T2 t) const
  {
    // m = (t mod R) * (-mod^(-1) mod R) mod R
    T m = (T)t * neg_r_inv; // (T)t performs t mod R; result m is also mod R
    // tmp = t + m*mod. tmp must be divisible by R.
    T2 tmp = t + (T2)m * mod;
    // res_intermediate = tmp / R. This result should be < 2*mod.
    T_mont res = (T)(tmp >> word_bits);
    // Ensure result is in [0, mod-1]
    res = res - (mod & -(res >= mod)); // if (res >= mod) res -= mod;
    return res;
  }

  T_mont add(T_mont a_mont, T_mont b_mont) const { return (a_mont + b_mont >= mod) ? (a_mont + b_mont - mod) : (a_mont + b_mont); }
  T_mont sub(T_mont a_mont, T_mont b_mont) const { return (a_mont >= b_mont) ? (a_mont - b_mont) : (a_mont + mod - b_mont); }
  T_mont mul(T_mont a_mont, T_mont b_mont) const { return reduce((T2)a_mont * b_mont); }
  /**
   * @brief 快速幂
   *
   * @param base_mont 底数（位于 Montgomery 数域）
   * @param exp 指数（普通整数）
   * @return T_mont 结果（位于 Montgomery 数域）
   */
  T_mont pow(T_mont base_mont, T exp) const
  {
    T_mont result_mont = from_T(1);
    while (exp > 0)
    {
      if (exp & 1)
        result_mont = mul(result_mont, base_mont);
      base_mont = mul(base_mont, base_mont);
      exp >>= 1;
    }
    return result_mont;
  }
  T_mont inv(T_mont x_mont) const { return pow(x_mont, mod - 2); }

private:
  T mod;       // 模数
  T r2;        // R^2 mod mod (where R = 2^mont_word_bits)
  T neg_r_inv; // -mod^(-1) mod R (where R = 2^mont_word_bits)
};
