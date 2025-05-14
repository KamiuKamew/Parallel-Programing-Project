#pragma once

#include <stdint.h>

using u32 = uint32_t;
using u64 = uint64_t;

using u32_mont = u32;
using u64_mont = u64;

/** 将 n 扩展为 2 的幂。 */
inline u32 expand_n(u32 n)
{
    u32 lg_n = 0;
    while ((1u << lg_n) < n)
        ++lg_n;
    return 1 << lg_n;
}

/** 将 a 的长度扩展为 2 的幂。 */
inline u32 *expand_a(u32 *a, u32 n, u32 n_expanded)
{
    // 虽然但是，这里应该有一个检查n是不是2的幂的逻辑
    // 话又说回来，这样得用别的库
    // 然而我又不想用别的库
    // 所以这里就先不检查了

    u32 *a_expanded = new u32[n_expanded];
    for (u32 i = 0; i < n; ++i)
        a_expanded[i] = a[i];
    for (u32 i = n; i < n_expanded; ++i)
        a_expanded[i] = 0;
    return a_expanded;
}

/**
 * @brief 就地对 a 做 bit-reverse 置换。
 *
 * 这个函数看起来不像是能 SIMD 优化的。
 * 而且这个函数应该在向量化之前用。
 *
 * @param a 输入序列（是不是 Montgomery 数域无所谓）
 * @param n 序列长度（扩展过，是2的幂）
 */
inline void bit_reverse_permute(u32 *a, u32 n)
{
    u32 lg_n = 0;
    while ((1u << lg_n) < n)
        ++lg_n;

    for (u32 i = 0; i < n; ++i)
    {
        u32 j = 0;
        for (u32 k = 0; k < lg_n; ++k)
            if (i & (1 << k))
                j |= (1 << (lg_n - 1 - k));
        if (i < j)
        {
            auto tmp = a[i];
            a[i] = a[j];
            a[j] = tmp;
        }
    }
}

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
inline void ntt_forward_mont(u32_mont *a_mont, u32 n, u32 p, u32_mont omega_mont)
{
    MontMod montMod(p);

    for (u32 mid = 1; mid < n; mid <<= 1)
    {
        u32_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1));
#pragma omp parallel for
        for (u32 j = 0; j < n; j += (mid << 1))
        {
            u32_mont w_mont = montMod.from_u32(1);
            for (u32 k = 0; k < mid; ++k, w_mont = montMod.mul(w_mont, Wn_mont))
            {
                u32_mont x_mont = a_mont[j + k];
                u32_mont y_mont = montMod.mul(w_mont, a_mont[j + k + mid]);
                a_mont[j + k] = montMod.add(x_mont, y_mont);
                a_mont[j + k + mid] = montMod.sub(x_mont, y_mont);
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
inline void ntt_inverse_mont(u32_mont *a_mont, u32 n, u32 p, u32_mont omega_mont)
{
    MontMod montMod(p);

    for (u32 mid = n >> 1; mid > 0; mid >>= 1)
    {
        u32_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1)); // Wn = ω⁻¹^((p-1)/(2*mid))
#pragma omp parallel for
        for (u32 j = 0; j < n; j += (mid << 1))
        {
            u32_mont w_mont = montMod.from_u32(1);
            for (u32 k = 0; k < mid; ++k, w_mont = montMod.mul(w_mont, Wn_mont))
            {
                u32_mont x_mont = a_mont[j + k];
                u32_mont y_mont = a_mont[j + k + mid];
                a_mont[j + k] = montMod.add(x_mont, y_mont);
                a_mont[j + k + mid] = montMod.mul(w_mont, montMod.sub(x_mont, y_mont));
            }
        }
    }

    u32_mont inv_n = montMod.inv(montMod.from_u32(n));
    for (u32 i = 0; i < n; ++i)
        a_mont[i] = montMod.mul(a_mont[i], inv_n);
}

#define OMEGA 3 // 998244353 的原根

/**
 * @brief 使用NTT优化的多项式乘法
 *
 * @param a 多项式系数
 * @param b 多项式系数
 * @param ab 结果
 * @param n 多项式长度
 * @param p 模数（质数）
 */
inline void poly_multiply_ntt(int *a, int *b, int *ab, int n, int p)
{
    MontMod montMod(p);

    u32 n_expanded = expand_n(2 * n - 1);
    u32 *a_expanded = expand_a((u32 *)a, n, n_expanded);
    u32 *b_expanded = expand_a((u32 *)b, n, n_expanded);

    bit_reverse_permute(a_expanded, n_expanded);
    bit_reverse_permute(b_expanded, n_expanded);

    u32_mont *a_mont = new u32_mont[n_expanded];
    u32_mont *b_mont = new u32_mont[n_expanded];
    u32_mont *ab_mont = new u32_mont[n_expanded];
    for (u32 i = 0; i < n_expanded; ++i)
        a_mont[i] = montMod.from_u32(a_expanded[i]);
    for (u32 i = 0; i < n_expanded; ++i)
        b_mont[i] = montMod.from_u32(b_expanded[i]);
    u32_mont omega_mont = montMod.from_u32(OMEGA);

    ntt_forward_mont(a_mont, n_expanded, p, omega_mont);
    ntt_forward_mont(b_mont, n_expanded, p, omega_mont);

    for (u32 i = 0; i < n_expanded; ++i)
        ab_mont[i] = montMod.mul(a_mont[i], b_mont[i]);

    ntt_inverse_mont(ab_mont, n_expanded, p, montMod.inv(omega_mont));

    for (u32 i = 0; i < n_expanded; ++i)
        ab[i] = montMod.to_u32(ab_mont[i]);

    bit_reverse_permute((u32 *)ab, n_expanded);
}