#pragma once

#include "general/op.h"

/**
 * @brief NTT 正变换：a(x) → A(ω)
 *
 * ntt_forward_mont 的 SIMD 版本。
 *
 * @param a_mont_simd 多项式系数（位于 Montgomery 数域），变换后表示频域系数（位于 Montgomery 数域）
 * @param n_expanded 多项式长度（普通整数）
 * @param p 模数（普通整数）
 * @param omega_mont 原根（位于 Montgomery 数域）
 */
inline void ntt_forward_mont_simd(u32x4_mont *a_mont_simd, u32 n_expanded, u32 p, u32_mont omega_mont)
{
}

/**
 * @brief NTT 逆变换：A(ω) → a_mont(x)
 *
 * ntt_inverse_dit_mont 的 SIMD 版本。
 *
 * @param a_mont_simd 频域系数，变换后表示多项式系数
 * @param n_expanded 多项式长度（普通整数）
 * @param p 模数（普通整数）
 * @param omega_mont 原根（位于 Montgomery 数域）
 */
inline void ntt_inverse_dit_mont_simd(u32x4_mont *a_mont_simd, u32 n_expanded, u32 p, u32_mont omega_mont)
{
}