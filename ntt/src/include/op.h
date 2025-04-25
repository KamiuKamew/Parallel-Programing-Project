#pragma once

#include "type.h"

// === 普通模运算 ===
inline u32 mod_add(u32 a, u32 b, u32 mod) { return (a + b) % mod; }
inline u32 mod_sub(u32 a, u32 b, u32 mod) { return (a - b + mod) % mod; }
inline u32 mod_mul(u32 a, u32 b, u32 mod) { return (1LL * a * b) % mod; }
inline u32 mod_pow(u32 base, u32 exp, u32 mod)
{
  u32 result = 1;
  while (exp > 0)
  {
    if (exp & 1)
    {
      result = mod_mul(result, base, mod);
    }
    base = mod_mul(base, base, mod);
    exp >>= 1;
  }
  return result;
}
inline u32 mod_inv(u32 x, u32 mod) { return mod_pow(x, mod - 2, mod); }
