#pragma once

#include "mod_base.h"

/**
 * @brief Barrett模运算实现
 *
 * 使用Barrett规约优化，在普通整数域中工作
 */
class ModBarrett : public ModBase
{
public:
    ModBarrett(u32 _mod, int iterations = 0) : ModBase(_mod)
    {
        barrett_k = 64;

        // 计算 r = floor(2^64 / mod)
        // 使用__uint128_t避免溢出
        __uint128_t power_of_2 = (__uint128_t)1 << barrett_k;
        u64 r = power_of_2 / _mod;

        // 迭代优化提高精度（可选）
        for (int i = 0; i < iterations; ++i)
        {
            __uint128_t numerator = power_of_2 - (__uint128_t)r * _mod;
            u64 correction = numerator / _mod;
            r += correction;
        }

        barrett_r = r;
    }

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
        u64 z = (u64)a * b;
        __uint128_t q_full = (__uint128_t)z * barrett_r;
        u64 q = q_full >> barrett_k;
        u32 res = (u32)(z - q * mod);
        if (res >= mod)
            res -= mod;
        if (res >= mod)
            res -= mod; // 再减一次以确保结果 < mod
        return res;
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

    // 数域转换实现（Barrett算法在普通域工作，无需转换）
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
        return "Barrett";
    }

private:
    u64 barrett_r; // 预计算的Barrett参数 r = floor(2^64 / mod)
    int barrett_k; // Barrett算法的k参数
};