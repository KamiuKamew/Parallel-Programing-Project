#pragma once

#include "type.h"

// === 普通模运算 ===
class Mod
{
public:
  Mod(u32 _mod) : mod(_mod) {}
  Mod(const Mod &) = delete;
  Mod &operator=(const Mod &) = delete;

  u32 add(u32 a, u32 b) const { return (a + b) % mod; }
  u32 sub(u32 a, u32 b) const { return (a >= b) ? (a - b) : (a + mod - b); }
  u32 mul(u32 a, u32 b) const { return (1LL * a * b) % mod; }
  u32 pow(u32 base, u32 exp) const
  {
    u32 result = 1;
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
  u32 inv(u32 x) const { return pow(x, mod - 2); }

private:
  u32 mod;
};

// === Montgomery 模运算 ===
class MontMod
{
public:
  MontMod(u32 _mod) : mod(_mod)
  {
    // 计算 r^2 mod mod，这里 r = 2^32
    u64 r2_temp = 1;
    for (int i = 0; i < 32; ++i)
      r2_temp = (r2_temp << 1) % mod;
    // 现在 r2_temp = 2^32 mod mod
    r2 = (u64)r2_temp * r2_temp % mod; // r^2 = (2^32)^2 mod mod

    // 计算 -mod^(-1) mod 2^32
    // 即满足 r*r_inv - mod*m_inv = 1 (mod r)，其中r=2^32
    u32 inv = 1;
    // 牛顿迭代法求模逆
    // 避免乘法溢出的方法是在32位计算中利用模2^32的特性
    for (int i = 0; i < 5; ++i)
      inv = (u32)((u64)inv * (2 - (u64)mod * inv));
    neg_r_inv = -inv;
  }

  MontMod(const MontMod &) = delete;
  MontMod &operator=(const MontMod &) = delete;

  u32_mont from_u32(u32 a) const { return reduce((u64)a * r2); }
  u32 to_u32(u32_mont a_mont) const { return reduce((u64)a_mont); }

  // Montgomery规约，计算 a * r^(-1) mod mod
  u32_mont reduce(u64 t) const
  {
    u32 m = (u32)t * neg_r_inv;
    u64 tmp = t + (u64)m * mod;
    u32_mont res = (u32)(tmp >> 32);
    res = res - (mod & -(res >= mod));
    return res;
  }

  u32_mont add(u32_mont a_mont, u32_mont b_mont) const { return (a_mont + b_mont >= mod) ? (a_mont + b_mont - mod) : (a_mont + b_mont); }
  u32_mont sub(u32_mont a_mont, u32_mont b_mont) const { return (a_mont >= b_mont) ? (a_mont - b_mont) : (a_mont + mod - b_mont); }
  u32_mont mul(u32_mont a_mont, u32_mont b_mont) const { return reduce((u64)a_mont * b_mont); }
  /**
   * @brief 快速幂
   *
   * @param base_mont 底数（位于 Montgomery 数域）
   * @param exp 指数（普通整数）
   * @return u32_mont 结果（位于 Montgomery 数域）
   */
  u32_mont pow(u32_mont base_mont, u32 exp) const
  {
    u32_mont result_mont = from_u32(1);
    while (exp > 0)
    {
      if (exp & 1)
        result_mont = mul(result_mont, base_mont);
      base_mont = mul(base_mont, base_mont);
      exp >>= 1;
    }
    return result_mont;
  }
  u32_mont inv(u32_mont x_mont) const { return pow(x_mont, mod - 2); }

private:
  u32 mod;       // 模数
  u32 r2;        // r^2 mod mod
  u32 neg_r_inv; // -r^(-1) mod 2^32
};
