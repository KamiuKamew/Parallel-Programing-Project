#pragma once

#include "general/type.h"

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