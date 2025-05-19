#pragma once

#include "../general/type.h"

inline u64 power_u64(u64 base, u64 exp, u64 mod)
{
    unsigned __int128 res_128 = 1;
    unsigned __int128 base_128 = base;

    base_128 %= mod;

    unsigned __int128 mod_128 = mod;

    while (exp > 0)
    {
        if (exp & 1)
            res_128 = (res_128 * base_128) % mod_128;
        base_128 = (base_128 * base_128) % mod_128;
        exp >>= 1;
    }
    return (u64)res_128;
}

// 通用CRT合并函数，支持任意多个模数
// residues: 各模数下的同一位置的余数数组
// mods: 各模数数组
// num_mods: 模数个数
// final_p: 最终取模
inline u32 crt_combine(const u32 *residues, const u32 *mods, int num_mods, u32 final_p)
{
    u64 P_prod = 1;
    for (int i = 0; i < num_mods; ++i)
        P_prod *= mods[i];

    u64 result = 0;
    for (int i = 0; i < num_mods; ++i)
    {
        u64 Mi = P_prod / mods[i];
        u64 ti = power_u64(Mi, mods[i] - 2, mods[i]);

        u64 term_contrib = residues[i]; // Implicitly u32 to u64
        term_contrib = static_cast<unsigned __int128>(term_contrib) * Mi % P_prod;
        term_contrib = static_cast<unsigned __int128>(term_contrib) * ti % P_prod;

        result = (static_cast<unsigned __int128>(result) + term_contrib) % P_prod;
    }

    return (u32)(result % final_p);
}
