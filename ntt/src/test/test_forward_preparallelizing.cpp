#include "../include/ntt.h"

#include <iostream>

int main()
{
    u32 a[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    u32 b[16] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};

    ntt_forward_mont(a, 16, 998244353, 3);
    ntt_forward_mont_before_simd(b, 16, 998244353, 3);

    for (int i = 0; i < 16; i++)
    {
        std::cout << a[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < 16; i++)
    {
        std::cout << b[i] << " ";
    }
    std::cout << std::endl;
}