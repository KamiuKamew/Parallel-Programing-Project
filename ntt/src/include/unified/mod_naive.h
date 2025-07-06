#pragma once

#include "mod_base.h"

/**
 * @brief 朴素模运算实现
 *
 * 使用直接的模运算，在普通整数域中工作
 */
class ModNaive : public ModBase
{
public:
    ModNaive(u32 _mod) : ModBase(_mod) {}

    // 基本模运算实现
    u32 add(u32 a, u32 b) const override
    {
        return (a + b) >= mod ? (a + b - mod) : (a + b);
    }

    u32 sub(u32 a, u32 b) const override
    {
        return (a >= b) ? (a - b) : (a + mod - b);
    }

    u32 mul(u32 a, u32 b) const override
    {
        return ((u64)a * b) % mod;
    }

    u32 pow(u32 base, u32 exp) const override
    {
        u32 result = 1;
        while (exp > 0)
        {
            if (exp & 1)
                result = mul(result, base);
            base = mul(base, base);
            exp >>= 1;
        }
        return result;
    }

    u32 inv(u32 x) const override
    {
        return pow(x, mod - 2);
    }

    // 数域转换实现（朴素算法在普通域工作，无需转换）
    u32 to_compute_domain(u32 a) const override
    {
        return a; // 无需转换
    }

    u32 from_compute_domain(u32 a) const override
    {
        return a; // 无需转换
    }

    const char *get_algorithm_name() const override
    {
        return "Naive";
    }
};