#pragma once

#include "../general/op.h"
#include "const.h"

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