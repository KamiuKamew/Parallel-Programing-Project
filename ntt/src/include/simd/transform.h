#pragma once

#include "../general/op.h"
#include "op.h"
#include "utils.h"
/**
 * @brief NTT 正变换：a(x) → A(ω)
 *
 * ntt_forward_mont 的 SIMD 版本。
 *
 * @param a_mont_simd 多项式系数（位于 Montgomery 数域），变换后表示频域系数（位于 Montgomery 数域）
 * @param n 多项式长度（普通整数）
 * @param p 模数（普通整数）
 * @param omega_mont 原根（位于 Montgomery 数域）
 */
inline void ntt_forward_mont_simd(u32x4_mont *a_mont_simd, u32 n, u32 p, u32_mont omega_mont)
{
    MontMod montMod(p);
    MontModNeon montModNeon(p);

    bool is_serial = false;
    u32_mont *a_mont = new u32_mont[n];

    for (u32 mid = 1; mid < n; mid <<= 1)
    {
        switch (mid)
        {
        case 1:
        {
            if (!is_serial)
            {
                from_simd(a_mont, a_mont_simd, n);
                is_serial = true;
            }
            for (u32 j = 0; j < n; j += (mid << 1))
            {
                u32_mont w_mont = montMod.from_u32(1);
                u32_mont x_mont = a_mont[j];
                u32_mont y_mont = montMod.mul(w_mont, a_mont[j + 1]);
                a_mont[j] = montMod.add(x_mont, y_mont);
                a_mont[j + 1] = montMod.sub(x_mont, y_mont);
            }
            break;
        }

        case 2:
        {
            if (!is_serial) // 实际上这不可能发生
            {
                from_simd(a_mont, a_mont_simd, n);
                is_serial = true;
            }
            u32_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1));
            for (u32 j = 0; j < n; j += (mid << 1))
            {
                u32_mont w_mont_0 = montMod.from_u32(1);
                u32_mont w_mont_1 = montMod.mul(w_mont_0, Wn_mont);
                u32_mont x_mont_0 = a_mont[j + 0];
                u32_mont x_mont_1 = a_mont[j + 1];
                u32_mont y_mont_0 = montMod.mul(w_mont_0, a_mont[j + 2]);
                u32_mont y_mont_1 = montMod.mul(w_mont_1, a_mont[j + 3]);
                a_mont[j + 0] = montMod.add(x_mont_0, y_mont_0);
                a_mont[j + 1] = montMod.add(x_mont_1, y_mont_1);
                a_mont[j + 2] = montMod.sub(x_mont_0, y_mont_0);
                a_mont[j + 3] = montMod.sub(x_mont_1, y_mont_1);
            }
            break;
        }

        default: // mid >= 4, parallelizable
        {
            if (is_serial)
            {
                to_simd(a_mont, a_mont_simd, n);
                is_serial = false;
            }
            u32_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1));
            for (u32 j = 0; j < n; j += (mid << 1))
            {
                u32_mont w_mont_0 = montMod.from_u32(1);
                u32_mont w_mont_1 = montMod.mul(w_mont_0, Wn_mont);
                u32_mont w_mont_2 = montMod.mul(w_mont_1, Wn_mont);
                u32_mont w_mont_3 = montMod.mul(w_mont_2, Wn_mont);
                u32_mont w_monts[4] = {w_mont_0, w_mont_1, w_mont_2, w_mont_3};
                u32x4_mont w_monts_simd = vld1q_u32(w_monts);

                u32_mont Wn_mont_4 = montMod.pow(Wn_mont, 4);
                u32x4_mont Wn_mont_4_simd = vdupq_n_u32(Wn_mont_4);

                for (u32 k = 0; k < mid; k += 4)
                {
                    // u32_mont x_mont_0 = a_mont[j + k + 0];
                    // u32_mont x_mont_1 = a_mont[j + k + 1];
                    // u32_mont x_mont_2 = a_mont[j + k + 2];
                    // u32_mont x_mont_3 = a_mont[j + k + 3];
                    u32x4_mont x_monts_simd = a_mont_simd[(j + k) / 4];

                    // u32_mont y_mont_0 = montMod.mul(w_mont_0, a_mont[j + k + mid + 0]);
                    // u32_mont y_mont_1 = montMod.mul(w_mont_1, a_mont[j + k + mid + 1]);
                    // u32_mont y_mont_2 = montMod.mul(w_mont_2, a_mont[j + k + mid + 2]);
                    // u32_mont y_mont_3 = montMod.mul(w_mont_3, a_mont[j + k + mid + 3]);
                    u32x4_mont y_monts_simd = montModNeon.mul(w_monts_simd, a_mont_simd[(j + k + mid) / 4]);

                    // a_mont[j + k + 0] = montMod.add(x_mont_0, y_mont_0);
                    // a_mont[j + k + 1] = montMod.add(x_mont_1, y_mont_1);
                    // a_mont[j + k + 2] = montMod.add(x_mont_2, y_mont_2);
                    // a_mont[j + k + 3] = montMod.add(x_mont_3, y_mont_3);
                    a_mont_simd[(j + k) / 4] = montModNeon.add(x_monts_simd, y_monts_simd);

                    // a_mont[j + k + mid + 0] = montMod.sub(x_mont_0, y_mont_0);
                    // a_mont[j + k + mid + 1] = montMod.sub(x_mont_1, y_mont_1);
                    // a_mont[j + k + mid + 2] = montMod.sub(x_mont_2, y_mont_2);
                    // a_mont[j + k + mid + 3] = montMod.sub(x_mont_3, y_mont_3);
                    a_mont_simd[(j + k + mid) / 4] = montModNeon.sub(x_monts_simd, y_monts_simd);

                    // w_mont_0 = montMod.mul(w_mont_0, Wn_mont_4);
                    // w_mont_1 = montMod.mul(w_mont_1, Wn_mont_4);
                    // w_mont_2 = montMod.mul(w_mont_2, Wn_mont_4);
                    // w_mont_3 = montMod.mul(w_mont_3, Wn_mont_4);
                    w_monts_simd = montModNeon.mul(w_monts_simd, Wn_mont_4_simd);
                }
            }
            break;
        }
        }
    }

    if (is_serial)
        to_simd(a_mont, a_mont_simd, n);

    delete[] a_mont;
}

/**
 * @brief NTT 逆变换：A(ω) → a_mont(x)
 *
 * ntt_inverse_dit_mont 的 SIMD 版本。
 *
 * @param a_mont_simd 频域系数，变换后表示多项式系数
 * @param n 多项式长度（普通整数）
 * @param p 模数（普通整数）
 * @param omega_mont 原根（位于 Montgomery 数域）
 */
inline void ntt_inverse_dit_mont_simd(u32x4_mont *a_mont_simd, u32 n, u32 p, u32_mont omega_mont)
{
}