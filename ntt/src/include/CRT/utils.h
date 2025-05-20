#pragma once

#include "constants.h"
#include <algorithm> // For std::swap, std::copy, std::fill
// #include <vector> // Not strictly needed for the functions being moved if using raw arrays

// 快速幂: 计算 (a^b) % p
inline u64 qpow(u64 a, u64 b, u64 p)
{
    u64 ans = 1;
    a %= p; // 确保底数在模p范围内
    while (b)
    {
        if (b & 1)
            ans = (1LL * ans * a) % p;
        a = (1LL * a * a) % p;
        b >>= 1;
    }
    return ans;
}

// 位逆序置换
inline void bit_reverse(u64 *a, int n)
{
    int j = 0;
    for (int i = 1; i < n; i++)
    {
        int bit = n >> 1;
        while (j & bit)
        {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j)
            std::swap(a[i], a[j]);
    }
}

// 扩展欧几里得算法: 计算 ax + by = gcd(a, b)
// 返回 gcd(a, b), x 和 y 通过引用返回
inline u64 extend_gcd(u64 a, u64 b, u64 &x, u64 &y)
{
    if (b == 0)
    {
        x = 1;
        y = 0;
        return a;
    }
    u64 x1, y1;
    u64 d = extend_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return d;
}

// 模逆元: 计算 a 在模 p 下的逆元
// 如果逆元不存在 (a 和 p 不互质)，行为未定义或返回错误值 (此处返回 -1)
inline u64 inv(u64 a, u64 p)
{
    u64 x, y;
    u64 d = extend_gcd(a, p, x, y);
    if (d != 1)
        return -1; // 逆元不存在
    return (x % p + p) % p;
}

// 特殊模乘: (a * b) % p，尝试避免溢出
inline u64 mulmod(u64 a, u64 b, u64 p)
{
    // unsigned long long x = 1Uu64 * a % p, y = 1Uu64 * b % p;
    // u64 result = x * y - (LL)((long double)x / p * y + 0.5L) * p;
    // return result < 0 ? result + p : result;
    // 使用更简单和标准的方法，前提是 a*b 不会溢出 long long
    a %= p;
    b %= p;
    unsigned __int128 res = (unsigned __int128)a * b;
    return res % p;
    // 如果 a*b 可能溢出 long long, 上面的 long double 方法或 __int128 (如果可用) 是备选方案
    // test_critical.cpp中的版本:
    // unsigned long long x = 1Uu64 * a % p, y = 1Uu64 * b % p;
    // u64 result = x * y - (LL)((long double)x / p * y + 0.5) * p;
    // return result < 0 ? result + p : result;
    // 这里使用 __int128 来获得更高的精度，如果编译器支持的话。
    // 如果不支持 __int128，则需要回退到 long double 方法或确保中间乘积在 u64 范围内。
    // 为了普适性，暂时采用 test_critical.cpp 中的版本，但需注意 long double 的精度问题。
    long long x = a % p;
    long long y = b % p;
    long long ret = x * y - (long long)((long double)x / p * y + 0.5L) * p;
    return ret < 0 ? ret + p : ret;
}

// 数论变换 (NTT)
inline void ntt(u64 *a, int n, u64 MOD, bool invert = false)
{
    bit_reverse(a, n);

    for (int len = 2; len <= n; len <<= 1)
    {
        u64 wn = qpow(ROOT, (MOD - 1) / len, MOD);
        if (invert)
            wn = inv(wn, MOD); // 使用模逆元

        for (int i = 0; i < n; i += len)
        {
            u64 w = 1;
            for (int j = 0; j < len / 2; j++)
            {
                u64 u = a[i + j];
                u64 v = (1LL * w * a[i + j + len / 2]) % MOD;
                a[i + j] = (u + v) % MOD;
                a[i + j + len / 2] = (u - v + MOD) % MOD;
                w = (1LL * w * wn) % MOD;
            }
        }
    }

    if (invert)
    {
        u64 inv_n = inv(n, MOD);
        for (int i = 0; i < n; i++)
            a[i] = (1LL * a[i] * inv_n) % MOD;
    }
}

// 中国剩余定理 (CRT): 合并两个同余方程组
// x === ab1[i] (mod p1)
// x === ab2[i] (mod p2)
// 结果保存在 ab1[i] 中，模数为 P = p1 * p2
inline void CRT(u64 *ab1, u64 *ab2, int n_coeffs, u64 p1, u64 p2)
{
    u64 P = 1LL * p1 * p2;
    u64 inv_p1_mod_p2 = inv(p1, p2);

    for (int i = 0; i < n_coeffs; i++)
    {
        // x = a1 + k*p1
        // a1 + k*p1 === a2 (mod p2)
        // k*p1 === a2 - a1 (mod p2)
        // k === (a2 - a1) * inv(p1, p2) (mod p2)
        u64 a1 = ab1[i];
        u64 a2 = ab2[i];
        u64 k = mulmod((a2 - a1 + p2) % p2, inv_p1_mod_p2, p2);
        ab1[i] = (a1 + mulmod(k, p1, P)) % P;
    }
    // 以下是 test_critical.cpp 中的 CRT 版本，逻辑略有不同，主要是如何计算合并系数。
    // 我将采用上面更常见的 CRT 构造方法。
    // 如果需要 test_critical.cpp 中的精确版本，我可以替换。
    /*
    u64 P_ = 1LL * p1 * p2;
    u64 P1_ = P_ / p1; // p2
    u64 P2_ = P_ / p2; // p1
    u64 inv_P1_mod_p1 = inv(P1_, p1);
    u64 inv_P2_mod_p2 = inv(P2_, p2);
    // u64 temp1 = mulmod(P1_, inv_P1_mod_p1, P_); // (p2 * inv(p2,p1)) % (p1*p2)
    // u64 temp2 = mulmod(P2_, inv_P2_mod_p2, P_); // (p1 * inv(p1,p2)) % (p1*p2)

    // Garner's algorithm style terms t1 = P2 * inv(P2,p1) = p1*inv(p1,p2) and t2 = P1*inv(P1,p2) = p2*inv(p2,p1)
    // x = v1*t1 + v2*t2 (mod P_)
    // test_critical.cpp's version: ab1[i] = (1LL * m1 * ab1[i] % P2_ * P1_ + (1LL * m2 * ab2[i] % P1_ * P2_)) % P_;
    // m1 = inv_P1_mod_p1, P2_ = p1, P1_ = p2
    // term1: inv(p2,p1) * ab1[i] % p1 * p2
    // m2 = inv_P2_mod_p2, P1_ = p2, P2_ = p1
    // term2: inv(p1,p2) * ab2[i] % p2 * p1
    // This looks like it's solving for x = C1*M1*ab1[i] + C2*M2*ab2[i] where M1=p2, M2=p1
    // And C1*M1 = 1 (mod p1), C2*M2 = 1 (mod p2)
    // C1 = inv(p2,p1), C2 = inv(p1,p2)
    // So, x = inv(p2,p1)*p2*ab1[i] + inv(p1,p2)*p1*ab2[i]
    // This is a valid Garner like construction as well.
    // The provided uncommented CRT in test_critical.cpp (line 376) is:
    // u64 temp1 = inv_P1 % P2_ * P1_; // inv_P1 = inv(P1_,p1) = inv(p2,p1); P2_ = p1; P1_ = p2 => inv(p2,p1) % p1 * p2
    // u64 temp2 = inv_P2 % P1_ * P2_; // inv_P2 = inv(P2_,p2) = inv(p1,p2); P1_ = p2; P2_ = p1 => inv(p1,p2) % p2 * p1
    // ab1[i] = (1LL * m1 * ab1[i] % P2_ * P1_ + (1LL * m2 * ab2[i] % P1_ * P2_)) % P_;
    // where m1 = temp1 / P1_ = inv(p2,p1) % p1 (this is inv(p2,p1) itself if result is < p1)
    // and m2 = temp2 / P2_ = inv(p1,p2) % p2
    // This formula is: ( (inv(p2,p1)*ab1[i]) % p1 * p2 + (inv(p1,p2)*ab2[i]) % p2 * p1 ) % P_
    // This is also a valid CRT formulation.
    // Sticking to the one from test_critical.cpp for closer porting:
    */
    u64 P_ = 1LL * p1 * p2;
    u64 P1_val = p2;                     // P / p1
    u64 P2_val = p1;                     // P / p2
    u64 inv_P1_mod_p1 = inv(P1_val, p1); // inv(p2, p1)
    u64 inv_P2_mod_p2 = inv(P2_val, p2); // inv(p1, p2)

    // Reconstruction terms based on test_critical.cpp's CRT logic for x = a (mod p1), x = b (mod p2)
    // x = (a * inv(p2,p1) * p2 + b * inv(p1,p2) * p1) % (p1*p2)
    // Let M1 = inv(p2,p1) * p2
    // Let M2 = inv(p1,p2) * p1
    u64 M1 = mulmod(inv_P1_mod_p1, P1_val, P_);
    u64 M2 = mulmod(inv_P2_mod_p2, P2_val, P_);

    for (int i = 0; i < n_coeffs; i++)
    {
        u64 term1 = mulmod(ab1[i], M1, P_);
        u64 term2 = mulmod(ab2[i], M2, P_);
        ab1[i] = (term1 + term2) % P_;
    }
}

// 使用NTT的多项式乘法 (单个模数)
inline void ntt_multiply(u64 *a, u64 *b, u64 *ab, int n_orig, u64 MOD)
{
    int size = 1;
    while (size < 2 * n_orig) // 需要补长到2的幂次，至少为2*n_orig - 1的下一个2的幂次
        size <<= 1;

    // 需要复制a和b到具有新长度的临时数组，或确保原始数组足够大且可以修改
    // 假设 a, b, ab 已经分配了足够的空间 (size)
    // 如果原始 a, b 的内容不应被修改或长度不足，则需要创建副本
    // test_critical.cpp 的 ntt_multiply 会补零。
    // 我们在此处假设调用者已处理或允许补零。
    // 如果 a 和 b 的实际内容长度是 n_orig，那么多出来的部分需要是0。
    // 这里我们显式补零，以符合 test_critical.cpp 中的做法。
    for (int i = n_orig; i < size; i++)
    {
        if (a + i < a + size)
            a[i] = 0; // 检查指针有效性以防万一，尽管通常size是目标大小
        if (b + i < b + size)
            b[i] = 0;
    }
    // 如果ab数组不是全零，也需要初始化。
    // for(int i=0; i<size; ++i) ab[i] = 0; // 通常ab是结果数组，不需要预先填0，除非累加
    // 但在NTT相乘后直接赋值，所以没问题

    ntt(a, size, MOD, false);
    ntt(b, size, MOD, false);

    for (int i = 0; i < size; i++)
        ab[i] = (1LL * a[i] * b[i]) % MOD;

    ntt(ab, size, MOD, true);
    // 结果多项式的长度是 2*n_orig - 1。ntt_multiply返回的ab数组长度是size。
    // 调用者需要知道实际需要的系数数量。
}
