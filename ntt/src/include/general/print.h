#pragma once

#include "type.h"

#include <iostream>

/**
 * @brief 调试打印标量数组内容
 *
 * @param a 数组指针
 * @param n 数组长度
 * @param prefix 打印前缀信息
 */
inline void debug_print_array(u32 *a, u32 n, const char *prefix = "")
{
    std::cout << prefix;
    for (u32 i = 0; i < n; ++i)
    {
        std::cout << a[i];
        if (i < n - 1)
            std::cout << ", ";
    }
    std::cout << std::endl;
}
