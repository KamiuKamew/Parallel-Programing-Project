#pragma once

#include "../general/utils.h"
#include "../general/op.h"
#include "transform.h"

#define OMEGA 3 // 998244353 的原根

template <typename T>
inline void poly_multiply_ntt_omp(T *a, T *b, T *ab, T n, T p)
{
    using T_mont = T;

    MontMod<T> montMod(p);

    T n_expanded = expand_n(2 * n - 1);
    T *a_expanded = expand_a((T *)a, n, n_expanded);
    T *b_expanded = expand_a((T *)b, n, n_expanded);

    bit_reverse_permute(a_expanded, n_expanded);
    bit_reverse_permute(b_expanded, n_expanded);

    T_mont *a_mont = new T_mont[n_expanded];
    T_mont *b_mont = new T_mont[n_expanded];
    T_mont *ab_mont = new T_mont[n_expanded];
    for (T i = 0; i < n_expanded; ++i)
        a_mont[i] = montMod.from_T(a_expanded[i]);
    for (T i = 0; i < n_expanded; ++i)
        b_mont[i] = montMod.from_T(b_expanded[i]);
    T_mont omega_mont = montMod.from_T(OMEGA);

    ntt_forward_mont_omp(a_mont, n_expanded, p, omega_mont);
    ntt_forward_mont_omp(b_mont, n_expanded, p, omega_mont);

    for (T i = 0; i < n_expanded; ++i)
        ab_mont[i] = montMod.mul(a_mont[i], b_mont[i]);

    ntt_inverse_mont_omp(ab_mont, n_expanded, p, montMod.inv(omega_mont));

    for (T i = 0; i < n_expanded; ++i)
        ab[i] = montMod.to_T(ab_mont[i]);

    bit_reverse_permute((T *)ab, n_expanded);
}