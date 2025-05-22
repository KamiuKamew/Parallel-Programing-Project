#pragma once

#include "../general/op.h"

template <typename T>
inline void ntt_forward_mont_omp(T *a_mont, T n, T p, T omega_mont)
{
    using T_mont = T;

    MontMod<T> montMod(p);

    for (T mid = 1; mid < n; mid <<= 1)
    {
        T_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1));
#pragma omp parallel for
        for (T j = 0; j < n; j += (mid << 1))
        {
            T_mont w_mont = montMod.from_T(1);
            for (T k = 0; k < mid; ++k, w_mont = montMod.mul(w_mont, Wn_mont))
            {
                T_mont x_mont = a_mont[j + k];
                T_mont y_mont = montMod.mul(w_mont, a_mont[j + k + mid]);
                a_mont[j + k] = montMod.add(x_mont, y_mont);
                a_mont[j + k + mid] = montMod.sub(x_mont, y_mont);
            }
        }
    }
}

template <typename T>
inline void ntt_inverse_mont_omp(T *a_mont, T n, T p, T omega_mont)
{
    using T_mont = T;

    MontMod<T> montMod(p);

    for (T mid = n >> 1; mid > 0; mid >>= 1)
    {
        T_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1)); // Wn = ω⁻¹^((p-1)/(2*mid))
#pragma omp parallel for
        for (T j = 0; j < n; j += (mid << 1))
        {
            T_mont w_mont = montMod.from_T(1);
            for (T k = 0; k < mid; ++k, w_mont = montMod.mul(w_mont, Wn_mont))
            {
                T_mont x_mont = a_mont[j + k];
                T_mont y_mont = a_mont[j + k + mid];
                a_mont[j + k] = montMod.add(x_mont, y_mont);
                a_mont[j + k + mid] = montMod.mul(w_mont, montMod.sub(x_mont, y_mont));
            }
        }
    }

    T_mont inv_n = montMod.inv(montMod.from_T(n));
    for (T i = 0; i < n; ++i)
        a_mont[i] = montMod.mul(a_mont[i], inv_n);
}