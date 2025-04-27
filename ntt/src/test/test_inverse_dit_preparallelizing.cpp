#include "../include/ntt.h"

#include <iostream>

int main()
{
    u32_mont a[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    u32_mont b[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};

    // 对使用相同的输入数据测试两个函数
    ntt_inverse_dit_mont(a, 16, 998244353, 3);
    ntt_inverse_mont_before_simd(b, 16, 998244353, 3);

    // 输出第一个函数的结果
    for (int i = 0; i < 16; i++)
    {
        std::cout << a[i] << " ";
    }
    std::cout << std::endl;

    // 输出第二个函数的结果
    for (int i = 0; i < 16; i++)
    {
        std::cout << b[i] << " ";
    }
    std::cout << std::endl;
}