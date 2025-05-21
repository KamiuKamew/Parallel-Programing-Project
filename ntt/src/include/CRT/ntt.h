/*
注：CRT部分仅支持T=u64, T2=u128，因此不进行模板化。
*/

#pragma once

#include "../general/utils.h"
#include "../general/op.h"
#include "../transform.h"

static const u64 CRT_MODS[] = {998244353, 1004535809, 469762049, 167772161};
static const u64 CRT_ROOTS[] = {3, 3, 3, 3};
static const u64 CRT_NUMS = sizeof(CRT_MODS) / sizeof(CRT_MODS[0]);

// 将 ab_crt 的 CRT 结果合并到 ab 中
inline void CRT_combine(u128 *ab, u64 **ab_crt, u64 n)
{

    u128 m = CRT_MODS[0];
    for (u64 i = 1; i < CRT_NUMS; ++i)
    {
        u64 CRT_MOD = CRT_MODS[i];

        Mod128 mod(CRT_MOD);

        u128 inv = mod.inv(m % CRT_MOD);
        for (u64 j = 0; j < n; j++)
        {
            u128 x = ab[j];
            u64 t = mod.sub(ab_crt[i][j], x % CRT_MOD);
            t = mod.mul(t, inv);

            x = x + m * t;
            ab[j] = x;
        }
        m *= CRT_MOD;
    }
}

inline void CRT_combine_2(u128 *ab, u64 **ab_crt, u64 n)
{
    for (u64 i = 0; i < n; ++i)
    {
        u128 x = ab_crt[0][i];
        u128 m = CRT_MODS[0];

        for (u64 j = 1; j < CRT_NUMS; ++j)
        {
            u64 CRT_MOD = CRT_MODS[j];
            Mod128 mod(CRT_MOD);

            u64 t = mod.sub(ab_crt[j][i], x % CRT_MOD);
            u64 inv = mod.inv(m % CRT_MOD);
            t = mod.mul(t, inv);

            x = x + m * t;
            m *= CRT_MOD;
        }

        ab[i] = x;
    }
}

/**
 * @brief 使用NTT和CRT优化的多项式乘法
 *
 * @param a 多项式系数 (assumed non-negative)
 * @param b 多项式系数 (assumed non-negative)
 * @param ab 结果多项式系数, modulo p
 * @param n 多项式长度 (number of coefficients, e.g., degree n-1)
 * @param p 最终结果的模数
 */
inline void poly_multiply_ntt_crt(u64 *a, u64 *b, u64 *ab, u64 n, u64 p)
{
    u64 n_expanded = expand_n(2 * n - 1);

    u64 **ab_crt = new u64 *[CRT_NUMS];
    u128 *ab_u128 = new u128[n_expanded];

    for (u64 i = 0; i < CRT_NUMS; i++)
    {
        ab_crt[i] = new u64[n_expanded]{};
        poly_multiply_ntt(a, b, ab_crt[i], n, CRT_MODS[i], CRT_ROOTS[i]);
    }

    for (u64 i = 0; i < n_expanded; ++i)
        ab_u128[i] = ab_crt[0][i];

    CRT_combine(ab_u128, ab_crt, n_expanded);
    // CRT_combine_garner(ab_u128, ab_crt, n_expanded);

    for (u64 i = 0; i < n_expanded; ++i)
        ab[i] = ab_u128[i] % p;

    delete[] ab_u128;
    for (u64 i = 0; i < CRT_NUMS; ++i)
        delete[] ab_crt[i];
    delete[] ab_crt;
}