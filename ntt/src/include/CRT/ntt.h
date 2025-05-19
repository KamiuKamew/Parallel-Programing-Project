#pragma once
#include "../general/utils.h"
#include "../transform.h"
#include "utils.h"
#include "constants.h"

/**
 * @brief 使用NTT和CRT优化的多项式乘法
 *
 * @param a 多项式系数 (assumed non-negative)
 * @param b 多项式系数 (assumed non-negative)
 * @param ab 结果多项式系数, modulo p
 * @param n 多项式长度 (number of coefficients, e.g., degree n-1)
 * @param p 最终结果的模数
 */
inline void poly_multiply_ntt_crt(int *a, int *b, int *ab, int n, int p)
{
    u32 n_expanded = expand_n(2 * n - 1);
    u32 *a_expanded = expand_a((u32 *)a, n, n_expanded);
    u32 *b_expanded = expand_a((u32 *)b, n, n_expanded);

    bit_reverse_permute(a_expanded, n_expanded);
    bit_reverse_permute(b_expanded, n_expanded);

    u32 **ab_bit_reversed = new u32 *[CRT_NUMS];
    for (int i = 0; i < CRT_NUMS; ++i)
    {

        MontMod montMod(CRT_MODS[i]);

        u32_mont *a_mont = new u32_mont[n_expanded];
        u32_mont *b_mont = new u32_mont[n_expanded];
        u32_mont *ab_mont = new u32_mont[n_expanded];
        for (u32 k = 0; k < n_expanded; ++k)
            a_mont[k] = montMod.from_u32(a_expanded[k] % CRT_MODS[i]);
        for (u32 k = 0; k < n_expanded; ++k)
            b_mont[k] = montMod.from_u32(b_expanded[k] % CRT_MODS[i]);
        u32_mont omega_mont = montMod.from_u32(CRT_ROOTS[i]);

        ntt_forward_mont(a_mont, n_expanded, CRT_MODS[i], omega_mont);
        ntt_forward_mont(b_mont, n_expanded, CRT_MODS[i], omega_mont);

        for (u32 k = 0; k < n_expanded; ++k)
            ab_mont[k] = montMod.mul(a_mont[k], b_mont[k]);

        ntt_inverse_mont(ab_mont, n_expanded, CRT_MODS[i], montMod.inv(omega_mont));

        ab_bit_reversed[i] = new u32[n_expanded];
        for (u32 k = 0; k < n_expanded; ++k)
            ab_bit_reversed[i][k] = montMod.to_u32(ab_mont[k]);

        delete[] a_mont;
        delete[] b_mont;
        delete[] ab_mont;
    }

    // CRT Reconstruction (通用多模数版本)
    u32 *ab_u32 = (u32 *)ab;               // operate on ab as u32 array
    u32 *crt_residues = new u32[CRT_NUMS]; // 临时存储每个模数下的同一位置余数
    for (u32 k = 0; k < n_expanded; ++k)
    {
        for (int i = 0; i < CRT_NUMS; ++i)
            crt_residues[i] = ab_bit_reversed[i][k];
        ab_u32[k] = crt_combine(crt_residues, CRT_MODS, CRT_NUMS, p);
    }
    delete[] crt_residues;

    bit_reverse_permute(ab_u32, n_expanded);

    // Cleanup
    delete[] a_expanded;
    delete[] b_expanded;
    for (int i = 0; i < CRT_NUMS; ++i)
        delete[] ab_bit_reversed[i];
    delete[] ab_bit_reversed;
}