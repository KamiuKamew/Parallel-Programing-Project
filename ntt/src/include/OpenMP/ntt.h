#pragma once

#include "../general/utils.h"
#include "../general/op.h"
#include "transform.h"

#define OMEGA 3 // 998244353 的原根

inline void poly_multiply_ntt_omp(int *a, int *b, int *ab, int n, int p)
{
    MontMod montMod(p);

    u32 n_expanded = expand_n(2 * n - 1);
    u32 *a_expanded = expand_a((u32 *)a, n, n_expanded);
    u32 *b_expanded = expand_a((u32 *)b, n, n_expanded);

    bit_reverse_permute(a_expanded, n_expanded);
    bit_reverse_permute(b_expanded, n_expanded);

    u32_mont *a_mont = new u32_mont[n_expanded];
    u32_mont *b_mont = new u32_mont[n_expanded];
    u32_mont *ab_mont = new u32_mont[n_expanded];
    for (u32 i = 0; i < n_expanded; ++i)
        a_mont[i] = montMod.from_u32(a_expanded[i]);
    for (u32 i = 0; i < n_expanded; ++i)
        b_mont[i] = montMod.from_u32(b_expanded[i]);
    u32_mont omega_mont = montMod.from_u32(OMEGA);

    ntt_forward_mont_omp(a_mont, n_expanded, p, omega_mont);
    ntt_forward_mont_omp(b_mont, n_expanded, p, omega_mont);

    for (u32 i = 0; i < n_expanded; ++i)
        ab_mont[i] = montMod.mul(a_mont[i], b_mont[i]);

    ntt_inverse_mont_omp(ab_mont, n_expanded, p, montMod.inv(omega_mont));

    for (u32 i = 0; i < n_expanded; ++i)
        ab[i] = montMod.to_u32(ab_mont[i]);

    bit_reverse_permute((u32 *)ab, n_expanded);
}