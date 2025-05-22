#pragma once

#include "../general/utils.h"
#include "utils.h"
#include "transform.h"

/**
 * @brief 使用NTT优化的多项式乘法
 *
 * 进一步实现了SIMD优化。
 *
 * @param a 多项式系数
 * @param b 多项式系数
 * @param ab 结果
 * @param n 多项式长度
 * @param p 模数（质数）
 */
inline void poly_multiply_ntt_simd(int *a, int *b, int *ab, int n, int p, int OMEGA = 3)
{
    MontModNeon montModNeon(p);
    MontMod<u32> montMod(p);

    u32 n_expanded = expand_n(2 * n - 1);
    u32 *a_expanded = expand_a<u32>((u32 *)a, (u32)n, n_expanded);
    u32 *b_expanded = expand_a<u32>((u32 *)b, (u32)n, n_expanded);

    bit_reverse_permute(a_expanded, n_expanded);
    bit_reverse_permute(b_expanded, n_expanded);

    u32x4 *a_simd = new u32x4[n_expanded / 4];
    u32x4 *b_simd = new u32x4[n_expanded / 4];
    u32x4 *ab_simd = new u32x4[n_expanded / 4];
    to_simd(a_expanded, a_simd, n_expanded);
    to_simd(b_expanded, b_simd, n_expanded);
    to_simd(new u32[n_expanded], ab_simd, n_expanded);
    u32 n_simd = n_expanded / 4;

    u32x4_mont *a_mont_simd = new u32x4_mont[n_simd];
    u32x4_mont *b_mont_simd = new u32x4_mont[n_simd];
    u32x4_mont *ab_mont_simd = new u32x4_mont[n_simd];
    for (u32 i = 0; i < n_simd; ++i)
        a_mont_simd[i] = montModNeon.from_u32x4(a_simd[i]);
    for (u32 i = 0; i < n_simd; ++i)
        b_mont_simd[i] = montModNeon.from_u32x4(b_simd[i]);
    u32_mont omega_mont = montMod.from_T(OMEGA);

    ntt_forward_mont_simd(a_mont_simd, n_expanded, p, omega_mont);
    ntt_forward_mont_simd(b_mont_simd, n_expanded, p, omega_mont);

    for (u32 i = 0; i < n_simd; ++i)
        ab_mont_simd[i] = montModNeon.mul(a_mont_simd[i], b_mont_simd[i]);

    ntt_inverse_mont_simd(ab_mont_simd, n_expanded, p, montMod.inv(omega_mont));

    for (u32 i = 0; i < n_simd; ++i)
        ab_simd[i] = montModNeon.to_u32x4(ab_mont_simd[i]); // 消除 mont

    u32 *ab_result = new u32[n_expanded];
    from_simd(ab_result, ab_simd, n_expanded); // 消除 simd

    bit_reverse_permute((u32 *)ab_result, n_expanded); // 消除 bit reversion

    for (u32 i = 0; i < n_expanded; ++i)
        ab[i] = ab_result[i]; // TODO: 这里可以跟上面的from simd合并。

    delete[] a_expanded;
    delete[] b_expanded;
    delete[] ab_result;
    delete[] a_mont_simd;
    delete[] b_mont_simd;
    delete[] ab_mont_simd;
}