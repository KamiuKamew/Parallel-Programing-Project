#pragma once

#include "general/type.h"

/** 将 a 转换成 SIMD 类型。 */
inline void to_simd(u32 *a, u32x4 *a_simd, u32 n_expanded)
{
    for (u32 i = 0; i < n_expanded / 4; ++i)
        a_simd[i] = vld1q_u32(&a[i * 4]);
}

/** 将 a_simd 转换成普通类型。 */
inline void from_simd(u32 *a, u32x4 *a_simd, u32 n_expanded)
{
    for (u32 i = 0; i < n_expanded / 4; ++i)
        vst1q_u32(&a[i * 4], a_simd[i]);
}

inline u32 get_lane(u32x4 v, int lane)
{
    switch (lane)
    {
    case 0:
        return vgetq_lane_u32(v, 0);
    case 1:
        return vgetq_lane_u32(v, 1);
    case 2:
        return vgetq_lane_u32(v, 2);
    case 3:
        return vgetq_lane_u32(v, 3);
    default:
        __builtin_unreachable();
    }
}

inline u32x4 set_lane(u32x4 v, u32 val, int lane)
{
    switch (lane)
    {
    case 0:
        return vsetq_lane_u32(val, v, 0);
    case 1:
        return vsetq_lane_u32(val, v, 1);
    case 2:
        return vsetq_lane_u32(val, v, 2);
    case 3:
        return vsetq_lane_u32(val, v, 3);
    default:
        __builtin_unreachable();
    }
}