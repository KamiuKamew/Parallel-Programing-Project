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
  T mul(T a, T b) const { return ((T2)a * (T2)b) % mod; }
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

class Mod128
{
public:
  Mod128(u128 _mod) : mod(_mod) {}
  Mod128(const Mod128 &) = delete;
  Mod128 &operator=(const Mod128 &) = delete;

  u128 add(u128 a, u128 b) const { return (a + b) % mod; }
  u128 sub(u128 a, u128 b) const { return (a >= b) ? (a - b) : (a + mod - b); }
  u128 mul(u128 a, u128 b) const { return (a * b) % mod; } // TODO：a * b 可能溢出
  u128 pow(u128 base, u128 exp) const
  {
    u128 result = 1;
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
  u128 inv(u128 x) const { return pow(x, mod - 2); }

private:
  u128 mod;
};

// === Montgomery 模运算 ===
template <typename T>
class MontMod
{
  using T_mont = T;
  using T2 = t_widen<T>;

public:
  MontMod(T _mod) : mod(_mod)
  {
    // 计算 r^2 mod mod，这里 r = 2^32
    T2 r2_temp = 1;
    for (int i = 0; i < 32; ++i)
      r2_temp = (r2_temp << 1) % mod;
    // 现在 r2_temp = 2^32 mod mod
    r2 = (T2)r2_temp * r2_temp % mod; // r^2 = (2^32)^2 mod mod

    // 计算 -mod^(-1) mod 2^32
    // 即满足 r*r_inv - mod*m_inv = 1 (mod r)，其中r=2^32
    T inv = 1;
    // 牛顿迭代法求模逆
    // 避免乘法溢出的方法是在32位计算中利用模2^32的特性
    for (int i = 0; i < 5; ++i)
      inv = (T)((T2)inv * (2 - (T2)mod * inv));
    neg_r_inv = -inv;
  }

  MontMod(const MontMod &) = delete;
  MontMod &operator=(const MontMod &) = delete;

  T_mont from_T(T a) const { return reduce((T2)a * r2); }
  T to_T(T_mont a_mont) const { return reduce((T2)a_mont); }

  // Montgomery规约，计算 a * r^(-1) mod mod
  T_mont reduce(T2 t) const
  {
    T m = (T)t * neg_r_inv;
    T2 tmp = t + (T2)m * mod;
    T_mont res = (T)(tmp >> 32);
    res = res - (mod & -(res >= mod));
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
  T r2;        // r^2 mod mod
  T neg_r_inv; // -r^(-1) mod 2^32
};
