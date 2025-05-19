#pragma once

#include "../general/op.h"

inline void ntt_forward_mont_omp(u32_mont *a_mont, u32 n, u32 p, u32_mont omega_mont)
{
    MontMod montMod(p);

    for (u32 mid = 1; mid < n; mid <<= 1)
    {
        u32_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1));
#pragma omp parallel for
        for (u32 j = 0; j < n; j += (mid << 1))
        {
            u32_mont w_mont = montMod.from_u32(1);
            for (u32 k = 0; k < mid; ++k, w_mont = montMod.mul(w_mont, Wn_mont))
            {
                u32_mont x_mont = a_mont[j + k];
                u32_mont y_mont = montMod.mul(w_mont, a_mont[j + k + mid]);
                a_mont[j + k] = montMod.add(x_mont, y_mont);
                a_mont[j + k + mid] = montMod.sub(x_mont, y_mont);
            }
        }
    }
}

inline void ntt_inverse_mont_omp(u32_mont *a_mont, u32 n, u32 p, u32_mont omega_mont)
{
    MontMod montMod(p);

    for (u32 mid = n >> 1; mid > 0; mid >>= 1)
    {
        u32_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1)); // Wn = ω⁻¹^((p-1)/(2*mid))
#pragma omp parallel for
        for (u32 j = 0; j < n; j += (mid << 1))
        {
            u32_mont w_mont = montMod.from_u32(1);
            for (u32 k = 0; k < mid; ++k, w_mont = montMod.mul(w_mont, Wn_mont))
            {
                u32_mont x_mont = a_mont[j + k];
                u32_mont y_mont = a_mont[j + k + mid];
                a_mont[j + k] = montMod.add(x_mont, y_mont);
                a_mont[j + k + mid] = montMod.mul(w_mont, montMod.sub(x_mont, y_mont));
            }
        }
    }

    u32_mont inv_n = montMod.inv(montMod.from_u32(n));
    for (u32 i = 0; i < n; ++i)
        a_mont[i] = montMod.mul(a_mont[i], inv_n);
}