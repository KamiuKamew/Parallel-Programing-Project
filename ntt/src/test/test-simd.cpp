#include <cstdio>
#include <cstdlib> // rand
#include <ctime>

#include "../include/ntt.h"

void test_montgomery_correctness()
{
    constexpr int N = 16;      // 测试N个数据（可以调大）
    const u32 mod = 998244353; // 测试模数（可以换成任意合法的）

    MontMod scalar(mod);
    MontModNeon simd(mod);

    // 准备随机数据
    u32 a_scalar[N], b_scalar[N];
    u32x4 a_simd[N / 4], b_simd[N / 4];

    srand((unsigned int)time(0));

    for (int i = 0; i < N; ++i)
    {
        a_scalar[i] = rand() % mod;
        b_scalar[i] = rand() % mod;
    }

    for (int i = 0; i < N; i += 4)
    {
        a_simd[i / 4] = vld1q_u32(&a_scalar[i]);
        b_simd[i / 4] = vld1q_u32(&b_scalar[i]);
    }

    // === 测试 from_u32 和 to_u32 ===
    for (int i = 0; i < N; i += 4)
    {
        u32x4 a_mont_simd = simd.from_u32x4(a_simd[i / 4]);
        u32x4 a_back_simd = simd.to_u32x4(a_mont_simd);

        for (int j = 0; j < 4; ++j)
        {
            u32 expected = scalar.to_u32(scalar.from_u32(a_scalar[i + j]));
            if (a_back_simd[j] != expected)
            {
                printf("from_u32/to_u32 mismatch at %d: got %u, expected %u\n", i + j, a_back_simd[j], expected);
                return;
            }
        }
    }

    // === 测试 add/sub/mul ===
    for (int i = 0; i < N; i += 4)
    {
        u32x4 a_mont_simd = simd.from_u32x4(a_simd[i / 4]);
        u32x4 b_mont_simd = simd.from_u32x4(b_simd[i / 4]);

        // scalar 版本
        u32 a_mont_scalar[4], b_mont_scalar[4];
        for (int j = 0; j < 4; ++j)
        {
            a_mont_scalar[j] = scalar.from_u32(a_scalar[i + j]);
            b_mont_scalar[j] = scalar.from_u32(b_scalar[i + j]);
        }

        // === 测试加法 ===
        u32x4 add_simd = simd.add(a_mont_simd, b_mont_simd);
        for (int j = 0; j < 4; ++j)
        {
            u32 expected = scalar.add(a_mont_scalar[j], b_mont_scalar[j]);
            if (add_simd[j] != expected)
            {
                printf("add mismatch at %d: got %u, expected %u\n", i + j, add_simd[j], expected);
                return;
            }
        }

        // === 测试减法 ===
        u32x4 sub_simd = simd.sub(a_mont_simd, b_mont_simd);
        for (int j = 0; j < 4; ++j)
        {
            u32 expected = scalar.sub(a_mont_scalar[j], b_mont_scalar[j]);
            if (sub_simd[j] != expected)
            {
                printf("sub mismatch at %d: got %u, expected %u\n", i + j, sub_simd[j], expected);
                return;
            }
        }

        // === 测试乘法 ===
        u32x4 mul_simd = simd.mul(a_mont_simd, b_mont_simd);
        for (int j = 0; j < 4; ++j)
        {
            u32 expected = scalar.mul(a_mont_scalar[j], b_mont_scalar[j]);
            if (mul_simd[j] != expected)
            {
                printf("mul mismatch at %d: got %u, expected %u\n", i + j, mul_simd[j], expected);
                return;
            }
        }
    }

    printf("All tests passed!\n");
}

int main()
{
    test_montgomery_correctness();
    return 0;
}
