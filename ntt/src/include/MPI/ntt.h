/*
注：CRT部分仅支持T=u64, T2=u128，因此不进行模板化。
*/

#pragma once

#include "../general/utils.h"
#include "../OpenMP_Barrett/ntt.h"

#include "../CRT/const.h"
#include "../CRT/crt.h"

/**
 * @brief 使用NTT和CRT优化的多项式乘法
 *
 * @param a 多项式系数 (assumed non-negative)
 * @param b 多项式系数 (assumed non-negative)
 * @param ab 结果多项式系数, modulo p
 * @param n 多项式长度 (number of coefficients, e.g., degree n-1)
 * @param p 最终结果的模数
 */
inline void poly_multiply_ntt_mpi(u64 *a, u64 *b, u64 *ab, u64 n, u64 p)
{
    u64 n_expanded = expand_n(2 * n - 1);
    u64 **ab_crt = new u64 *[CRT_NUMS];
    u128 *ab_u128 = new u128[n_expanded];

    for (u64 i = 0; i < CRT_NUMS; i++)
    {
        ab_crt[i] = new u64[n_expanded]{};

        // 创建临时数组，对输入系数进行模数归约
        u64 *a_mod = new u64[n];
        u64 *b_mod = new u64[n];
        for (u64 j = 0; j < n; j++)
        {
            a_mod[j] = a[j] % CRT_MODS[i];
            b_mod[j] = b[j] % CRT_MODS[i];
        }

        /*
        ## 问题原因与解决方案

        ### **问题根本原因：**
        `poly_multiply_ntt_mpi` 与 `poly_multiply_ntt_crt` 的关键区别在于**输入系数的处理方式**：

        1. **`poly_multiply_ntt_crt`：** 使用模板版本的 `poly_multiply_ntt`，该函数内部会自动对输入系数进行模数归约
        2. **`poly_multiply_ntt_mpi`：** 直接调用 `poly_multiply_ntt_omp_Barrett`，但没有对输入系数进行预处理

        ### **具体问题：**
        当输入多项式系数超过CRT模数时（如大模数测试用例中的 `992009819965388`），直接传递给 `poly_multiply_ntt_omp_Barrett` 会导致：
        - 输入值远超过工作模数（如 `998244353`）
        - Barrett约简和Montgomery算法无法正确处理这些超大输入
        - 最终导致计算结果错误

        ### **解决方案：**
        在调用 `poly_multiply_ntt_omp_Barrett` 之前，对输入系数进行模数归约：

        ```cpp
        // 创建临时数组，对输入系数进行模数归约
        u64 *a_mod = new u64[n];
        u64 *b_mod = new u64[n];
        for (u64 j = 0; j < n; j++) {
            a_mod[j] = a[j] % CRT_MODS[i];
            b_mod[j] = b[j] % CRT_MODS[i];
        }

        poly_multiply_ntt_omp_Barrett(a_mod, b_mod, ab_crt[i], n, CRT_MODS[i], CRT_ROOTS[i]);
        ```

        ### **验证结果：**
        - ✅ 所有5个测试用例都通过
        - ✅ 包括最困难的 `n=131072, p=1337006139375617` 大模数测试
        - ✅ 性能保持在合理范围内（约500微秒）

        这个修复确保了 `poly_multiply_ntt_mpi` 能够正确处理任意大小的输入系数，使其行为与 `poly_multiply_ntt_crt` 完全一致。
        */

        poly_multiply_ntt_omp_Barrett(a_mod, b_mod, ab_crt[i], n, CRT_MODS[i], CRT_ROOTS[i]);

        delete[] a_mod;
        delete[] b_mod;
    }

    for (u64 i = 0; i < n_expanded; ++i)
        ab_u128[i] = ab_crt[0][i];

    CRT_combine(ab_u128, ab_crt, n_expanded);
    // CRT_combine_garner(ab_u128, ab_crt, n_expanded);

    for (u64 i = 0; i < n_expanded; ++i)
        ab[i] = ab_u128[i] % p;

    delete[] ab_u128;
    for (u64 i = 0; i < CRT_NUMS; ++i)
        delete[] ab_crt[i];
    delete[] ab_crt;
}