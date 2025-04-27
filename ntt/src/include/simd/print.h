#pragma once

#include "type.h"
#include "../general/type.h"

#include <iostream>

/**
 * @brief 调试打印SIMD向量数组内容
 *
 * @param a_simd SIMD向量数组指针
 * @param n_simd 向量数组长度
 * @param prefix 打印前缀信息
 */
inline void debug_print_simd_array(u32x4 *a_simd, u32 n_simd, const char *prefix = "")
{
    std::cout << prefix;
    for (u32 i = 0; i < n_simd; ++i)
    {
        std::cout << "[";
        std::cout << vgetq_lane_u32(a_simd[i], 0) << ", "
                  << vgetq_lane_u32(a_simd[i], 1) << ", "
                  << vgetq_lane_u32(a_simd[i], 2) << ", "
                  << vgetq_lane_u32(a_simd[i], 3);
        std::cout << "]";
        if (i < n_simd - 1)
            std::cout << ", ";
    }
    std::cout << std::endl;
}