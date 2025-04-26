#pragma once

#include "../general/type.h"
#include "type.h"

// === Montgomery 模运算（NEON） ===
class MontModNeon
{
public:
    MontModNeon(u32 _mod) : mod(_mod)
    {
        // 计算 r^2 mod mod，这里 r = 2^32
        u64 r2_temp = 1;
        for (int i = 0; i < 32; ++i)
            r2_temp = (r2_temp << 1) % mod;
        r2 = (u64)r2_temp * r2_temp % mod;

        u32 inv = 1;
        for (int i = 0; i < 5; ++i)
            inv = (u32)((u64)inv * (2 - (u64)mod * inv));
        neg_r_inv = -inv;

        mod_vec = vdupq_n_u32(mod);
        neg_r_inv_vec = vdupq_n_u32(neg_r_inv);
    }

    MontModNeon(const MontModNeon &) = delete;
    MontModNeon &operator=(const MontModNeon &) = delete;

    u32x4_mont from_u32x4(u32x4 a) const
    {
        u64x2 t0 = vmull_u32(vget_low_u32(a), vdup_n_u32(r2));
        u64x2 t1 = vmull_u32(vget_high_u32(a), vdup_n_u32(r2));
        return reduce_pair(t0, t1);
    }
    u32x4 to_u32x4(u32x4_mont a_mont) const { return reduce(a_mont); }

    // SIMD版 Montgomery规约
    u32x4_mont reduce(u32x4 t_lo) const
    {
        // t_lo 是低32位，高32位补0，提升成64位
        u64x2 t0 = vmovl_u32(vget_low_u32(t_lo));
        u64x2 t1 = vmovl_u32(vget_high_u32(t_lo));

        return reduce_pair(t0, t1);
    }

    u32x4_mont reduce_pair(u64x2 t0, u64x2 t1) const
    {
        // m = (t mod 2^32) * neg_r_inv mod 2^32
        u32x2 m0 = vmul_u32(vmovn_u64(t0), vget_low_u32(neg_r_inv_vec));
        u32x2 m1 = vmul_u32(vmovn_u64(t1), vget_high_u32(neg_r_inv_vec));

        // t + m * mod
        u64x2 t0_new = vmlal_u32(t0, m0, vget_low_u32(mod_vec));
        u64x2 t1_new = vmlal_u32(t1, m1, vget_high_u32(mod_vec));

        // (t + m * mod) >> 32
        u32x2 res0 = vshrn_n_u64(t0_new, 32);
        u32x2 res1 = vshrn_n_u64(t1_new, 32);

        u32x4 res = vcombine_u32(res0, res1);

        // res = res - (mod & -(res >= mod))
        uint32x4_t mask = vcgeq_u32(res, mod_vec); // res >= mod
        uint32x4_t mod_masked = vandq_u32(mod_vec, mask);
        res = vsubq_u32(res, mod_masked);

        return res;
    }

    u32x4_mont add(u32x4_mont a, u32x4_mont b) const
    {
        u32x4 res = vaddq_u32(a, b);
        uint32x4_t mask = vcgeq_u32(res, mod_vec);
        uint32x4_t mod_masked = vandq_u32(mod_vec, mask);
        return vsubq_u32(res, mod_masked);
    }

    u32x4_mont sub(u32x4_mont a, u32x4_mont b) const
    {
        uint32x4_t mask = vcgeq_u32(a, b);
        u32x4 res1 = vsubq_u32(a, b);
        u32x4 res2 = vsubq_u32(vaddq_u32(a, mod_vec), b);
        return vbslq_u32(mask, res1, res2);
    }

    u32x4_mont mul(u32x4_mont a, u32x4_mont b) const
    {
        u64x2 prod0 = vmull_u32(vget_low_u32(a), vget_low_u32(b));
        u64x2 prod1 = vmull_u32(vget_high_u32(a), vget_high_u32(b));
        return reduce_pair(prod0, prod1);
    }

    u32x4_mont pow(u32x4_mont base_mont, u32 exp) const
    {
        u32x4 result_mont = from_u32x4(vdupq_n_u32(1));
        while (exp > 0)
        {
            if (exp & 1)
                result_mont = mul(result_mont, base_mont);
            base_mont = mul(base_mont, base_mont);
            exp >>= 1;
        }
        return result_mont;
    }

    u32x4_mont inv(u32x4_mont x_mont) const { return pow(x_mont, mod - 2); }

private:
    u32 mod;
    u32 r2;
    u32 neg_r_inv;

    u32x4 mod_vec;       // 向量化的模数
    u32x4 neg_r_inv_vec; // 向量化的 -r^(-1) mod 2^32
};