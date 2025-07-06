#pragma once

#include "mod_base.h"

/**
 * @brief Montgomery模运算实现
 *
 * 使用Montgomery规约优化，在Montgomery域中工作
 */
class ModMontgomery : public ModBase
{
public:
    ModMontgomery(u32 _mod) : ModBase(_mod)
    {
        // 计算 R^2 mod mod，这里 R = 2^32
        u64 r_val_mod_n = 1;
        for (int i = 0; i < 32; ++i)
            r_val_mod_n = (r_val_mod_n << 1) % mod;
        r2 = (u64)r_val_mod_n * r_val_mod_n % mod;

        // 计算 -mod^(-1) mod R (R = 2^32)
        u32 inv = 1;
        for (int i = 0; i < 5; ++i) // 5次迭代足够32位
            inv = (u32)((u64)inv * (2 - (u64)mod * inv));
        neg_r_inv = -inv;
    }

    // 基本模运算实现（在Montgomery域中）
    u32 add(u32 a_mont, u32 b_mont) const override
    {
        return (a_mont + b_mont >= mod) ? (a_mont + b_mont - mod) : (a_mont + b_mont);
    }

    u32 sub(u32 a_mont, u32 b_mont) const override
    {
        return (a_mont >= b_mont) ? (a_mont - b_mont) : (a_mont + mod - b_mont);
    }

    u32 mul(u32 a_mont, u32 b_mont) const override
    {
        return mont_reduce((u64)a_mont * b_mont);
    }

    u32 pow(u32 base_mont, u32 exp) const override
    {
        u32 result_mont = to_compute_domain(1); // 1转换到Montgomery域
        while (exp > 0)
        {
            if (exp & 1)
                result_mont = mul(result_mont, base_mont);
            base_mont = mul(base_mont, base_mont);
            exp >>= 1;
        }
        return result_mont;
    }

    u32 inv(u32 x_mont) const override
    {
        return pow(x_mont, mod - 2);
    }

    // 数域转换实现
    u32 to_compute_domain(u32 a) const override
    {
        return mont_reduce((u64)a * r2); // 普通域 -> Montgomery域
    }

    u32 from_compute_domain(u32 a_mont) const override
    {
        return mont_reduce((u64)a_mont); // Montgomery域 -> 普通域
    }

    const char *get_algorithm_name() const override
    {
        return "Montgomery";
    }

private:
    u32 r2;        // R^2 mod mod (where R = 2^32)
    u32 neg_r_inv; // -mod^(-1) mod R (where R = 2^32)

    // Montgomery规约，计算 t * R^(-1) mod mod, where R = 2^32
    u32 mont_reduce(u64 t) const
    {
        u32 m = (u32)t * neg_r_inv;          // m = (t mod R) * (-mod^(-1) mod R) mod R
        u64 tmp = t + (u64)m * mod;          // tmp = t + m*mod
        u32 res = (u32)(tmp >> 32);          // res = tmp / R
        return res >= mod ? res - mod : res; // 确保结果在[0, mod)范围内
    }
};