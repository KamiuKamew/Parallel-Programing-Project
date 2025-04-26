#pragma once

#include "type.h"

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