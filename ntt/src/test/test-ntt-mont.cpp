#include "../include/ntt.h"

#include <iostream>

int main()
{
    u32 p = 998244353;
    MontMod montMod(p);

    u32 a[] = {1, 2, 3, 4};
    u32_mont a_mont[4];
    for (int i = 0; i < 4; ++i)
        a_mont[i] = montMod.from_u32(a[i]);

    // ntt forward
    std::cout << "ntt forward" << std::endl;

    ntt_forward(a, 4, p, 3);
    ntt_forward_mont(a_mont, 4, p, montMod.from_u32(3));

    u32 a_inv[4];
    for (int i = 0; i < 4; ++i)
        a_inv[i] = montMod.to_u32(a_mont[i]);

    for (int i = 0; i < 4; ++i)
        std::cout << a[i] << " ";
    std::cout << std::endl;
    for (int i = 0; i < 4; ++i)
        std::cout << a_inv[i] << " ";
    std::cout << std::endl;

    // ntt inverse
    std::cout << "ntt inverse" << std::endl;

    ntt_inverse(a, 4, p, 3);
    ntt_inverse_mont(a_mont, 4, p, montMod.from_u32(3));

    for (int i = 0; i < 4; ++i)
        a_inv[i] = montMod.to_u32(a_mont[i]);

    for (int i = 0; i < 4; ++i)
        std::cout << a[i] << " ";
    std::cout << std::endl;
    for (int i = 0; i < 4; ++i)
        std::cout << a_inv[i] << " ";
    std::cout << std::endl;
}
