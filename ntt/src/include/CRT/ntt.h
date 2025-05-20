#pragma once
#include "../general/utils.h"
#include "../transform.h"
#include "utils.h"
#include "constants.h"
#include <algorithm> // For std::copy, std::fill
#include <vector>    // For std::vector if used, or just for completeness if new/delete used

/**
 * @brief 使用NTT和CRT优化的多项式乘法
 *
 * @param a 多项式系数 (assumed non-negative, type LL*)
 * @param b 多项式系数 (assumed non-negative, type LL*)
 * @param ab 结果多项式系数, modulo p (type LL*)
 * @param n 多项式长度 (number of coefficients, e.g., degree n-1)
 * @param p 最终结果的模数 (type LL)
 */
inline void poly_multiply_ntt_crt(int *a, int *b, int *ab, int n, int p_final_mod)
{
    // 计算长度，用于创建复制变量
    int size = 1;
    while (size < 2 * n) // 结果多项式的系数个数最多为 2n-1, NTT长度需要是 >= 2n-1 的2的幂
        size <<= 1;

    // u64 *a_copy[NUM_NTT_PRIMES]; // 使用常量 NUM_NTT_PRIMES
    // u64 *b_copy[NUM_NTT_PRIMES];
    // u64 *ab_copy[NUM_NTT_PRIMES];
    // 使用 std::vector<LL*> 或者直接数组，但注意数组大小需为常量
    // test_critical.cpp uses 4 primes hardcoded
    const int num_primes_to_use = 4; // Corresponds to NTT_PRIMES in constants.h
    if (NUM_NTT_PRIMES < num_primes_to_use)
    {
        // Handle error: not enough primes defined
        // For now, assume NUM_NTT_PRIMES >= num_primes_to_use
    }

    std::vector<u64 *> a_copy(num_primes_to_use);
    std::vector<u64 *> b_copy(num_primes_to_use);
    std::vector<u64 *> ab_copy(num_primes_to_use);

    a_copy[0] = (u64 *)a; // 使用原始数组作为第一个副本，避免额外分配和复制
    b_copy[0] = (u64 *)b;
    ab_copy[0] = (u64 *)ab; // 结果将直接写入 ab

    for (int i = 1; i < num_primes_to_use; i++)
    {
        a_copy[i] = new u64[size];
        b_copy[i] = new u64[size];
        ab_copy[i] = new u64[size];
    }

    // 复制原始系数到其他副本 (从 a_copy[0] 和 b_copy[0] 复制)
    for (int i = 1; i < num_primes_to_use; i++)
    {
        std::copy(a_copy[0], a_copy[0] + n, a_copy[i]); // 只复制n个系数
        std::fill(a_copy[i] + n, a_copy[i] + size, 0);  // 其余补零
        std::copy(b_copy[0], b_copy[0] + n, b_copy[i]); // 只复制n个系数
        std::fill(b_copy[i] + n, b_copy[i] + size, 0);  // 其余补零
        // std::fill(ab_copy[i], ab_copy[i] + size, 0); // ntt_multiply内部会处理ab的初始化
    }
    // 确保 a_copy[0] 和 b_copy[0] 也补零 (如果n < size)
    std::fill(a_copy[0] + n, a_copy[0] + size, 0);
    std::fill(b_copy[0] + n, b_copy[0] + size, 0);

    // 对每个模数执行NTT乘法
    for (int i = 0; i < num_primes_to_use; i++)
    {
        ntt_multiply(a_copy[i], b_copy[i], ab_copy[i], n, NTT_PRIMES[i]);
    }

    // 使用CRT合并结果
    // 第一轮合并: (ab_copy[0], ab_copy[1]) -> ab_copy[0], (ab_copy[2], ab_copy[3]) -> ab_copy[2]
    for (int i = 0; i < num_primes_to_use / 2; i++)
    {
        CRT(ab_copy[2 * i], ab_copy[2 * i + 1], size, NTT_PRIMES[2 * i], NTT_PRIMES[2 * i + 1]);
    }

    // 第二轮合并 (Garner's algorithm like step for final modulus)
    // 当前 ab_copy[0] 中的结果模 P1 = NTT_PRIMES[0] * NTT_PRIMES[1]
    // 当前 ab_copy[2] 中的结果模 P2 = NTT_PRIMES[2] * NTT_PRIMES[3]
    u64 P1_crt = 1LL * NTT_PRIMES[0] * NTT_PRIMES[1];
    u64 P2_crt = 1LL * NTT_PRIMES[2] * NTT_PRIMES[3];

    // 我们需要解系统:
    // x === ab_copy[0][i] (mod P1_crt)
    // x === ab_copy[2][i] (mod P2_crt)
    // 并最终得到 x (mod p_final_mod)
    // x = ab_copy[0][i] + k * P1_crt
    // ab_copy[0][i] + k * P1_crt === ab_copy[2][i] (mod P2_crt)
    // k * P1_crt === ab_copy[2][i] - ab_copy[0][i] (mod P2_crt)
    // k === (ab_copy[2][i] - ab_copy[0][i]) * inv(P1_crt, P2_crt) (mod P2_crt)
    u64 inv_P1_mod_P2 = inv(P1_crt, P2_crt);
    u64 M_total_crt = 1LL * P1_crt * P2_crt; // 这是能表示的最大模数

    for (int i = 0; i < size; i++)
    {
        u64 r1 = ab_copy[0][i]; // mod P1_crt
        u64 r2 = ab_copy[2][i]; // mod P2_crt

        u64 k = mulmod((r2 - r1 + P2_crt) % P2_crt, inv_P1_mod_P2, P2_crt);
        u64 combined_res = (r1 + mulmod(k, P1_crt, M_total_crt)) % M_total_crt;
        ab[i] = combined_res % p_final_mod; // 最后取模到目标模数 p_final_mod
    }

    // 清理动态分配的内存
    for (int i = 1; i < num_primes_to_use; i++)
    {
        delete[] a_copy[i];
        delete[] b_copy[i];
        delete[] ab_copy[i];
    }
}