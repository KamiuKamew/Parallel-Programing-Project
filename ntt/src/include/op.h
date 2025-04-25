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
