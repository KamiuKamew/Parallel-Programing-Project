#pragma once

#include "../general/utils.h"
#include "../general/op.h"
#include "../transform.h"

static const int CRT_MODS[] = {998244353, 1004535809, 469762049};
static const int CRT_ROOTS[] = {3, 3, 3};
static const int CRT_NUMS = sizeof(CRT_MODS) / sizeof(CRT_MODS[0]);

// 将 ab_crt 的 CRT 结果合并到 ab 中
inline void CRT_combine(u128 *ab, int *ab_crt, int n, u128 p_ab, u128 p_crt)
{
    Mod128 mod_crt(p_crt);
    Mod128 mod_M(p_ab * p_crt);

    u128 inv = mod_crt.inv(p_ab);
    for (int i = 0; i < n; i++)
        ab[i] = mod_M.add(ab[i], mod_M.mul(p_ab, mod_crt.mul(mod_crt.sub(ab_crt[i], ab[i] % p_crt), inv)));
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
inline void poly_multiply_ntt_crt(int *a, int *b, int *ab, int n, int p)
{
    u32 n_expanded = expand_n(2 * n - 1);

    int **ab_crt = new int *[CRT_NUMS];
    u128 *ab_u128 = new u128[n_expanded];
    for (int i = 0; i < CRT_NUMS; i++)
    {
        ab_crt[i] = new int[n_expanded]{};
        poly_multiply_ntt(a, b, ab_crt[i], n, CRT_MODS[i], CRT_ROOTS[i]);
    }

    for (int i = 0; i < n_expanded; ++i)
        ab_u128[i] = ab_crt[0][i];

    u128 m = CRT_MODS[0];
    for (int i = 1; i < CRT_NUMS; ++i)
    {
        CRT_combine(ab_u128, ab_crt[i], n_expanded, m, CRT_MODS[i]);
        m *= CRT_MODS[i];
    }

    for (int i = 0; i < n_expanded; ++i)
        ab[i] = ab_u128[i] % p;

    delete[] ab_u128;
    for (int i = 0; i < CRT_NUMS; ++i)
        delete[] ab_crt[i];
    delete[] ab_crt;
}