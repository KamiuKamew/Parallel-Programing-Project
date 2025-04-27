#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <random>
#include <algorithm>
#include <vector>

#include "../include/general/type.h"
#include "../include/general/utils.h"
#include "../include/general/op.h"
#include "../include/transform.h"
#include "../include/simd/transform.h"

// 随机生成测试数据
void generate_random_data(u32 *data, u32 n, u32 p)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<u32> distrib(0, p - 1);

    for (u32 i = 0; i < n; ++i)
    {
        data[i] = distrib(gen);
    }
}

// 比较两个数组是否相等
bool compare_arrays(u32 *a, u32 *b, u32 n)
{
    for (u32 i = 0; i < n; ++i)
    {
        if (a[i] != b[i])
        {
            return false;
        }
    }
    return true;
}

// 测试ntt_inverse_dit_mont和ntt_inverse_dit_mont_simd的一致性
void test_ntt_inverse_consistency(u32 n, u32 p, u32_mont omega_mont)
{
    printf("测试 n=%u, p=%u 的NTT逆变换的一致性\n", n, p);

    // 确保n是2的幂
    u32 n_expanded = expand_n(n);

    // 生成随机数据
    u32 *a = new u32[n_expanded];
    generate_random_data(a, n_expanded, p);

    // 创建两个Montgomery域的副本进行测试
    MontMod montMod(p);
    u32_mont *a_mont1 = new u32_mont[n_expanded];
    u32_mont *a_mont2 = new u32_mont[n_expanded];

    // 转换为Montgomery域
    for (u32 i = 0; i < n_expanded; ++i)
    {
        a_mont1[i] = montMod.from_u32(a[i]);
        a_mont2[i] = montMod.from_u32(a[i]);
    }

    // 使用常规的NTT逆变换
    ntt_inverse_dit_mont(a_mont1, n_expanded, p, omega_mont);

    // 创建SIMD版本的数据
    u32_mont *a_mont2_copy = new u32_mont[n_expanded];
    std::copy(a_mont2, a_mont2 + n_expanded, a_mont2_copy);

    u32x4_mont *a_mont_simd = new u32x4_mont[n_expanded / 4];
    for (u32 i = 0; i < n_expanded / 4; ++i)
    {
        a_mont_simd[i] = vdupq_n_u32(0);
        for (u32 j = 0; j < 4; ++j)
        {
            a_mont_simd[i] = set_lane(a_mont_simd[i], a_mont2_copy[i * 4 + j], j);
        }
    }

    // 使用SIMD版本的NTT逆变换
    ntt_inverse_dit_mont_simd(a_mont_simd, n_expanded, p, omega_mont);

    // 将SIMD结果转换回普通数组
    u32_mont *a_mont_simd_result = new u32_mont[n_expanded];
    for (u32 i = 0; i < n_expanded / 4; ++i)
    {
        for (u32 j = 0; j < 4; ++j)
        {
            a_mont_simd_result[i * 4 + j] = get_lane(a_mont_simd[i], j);
        }
    }

    // 转换回普通域以便打印比较
    u32 *result1 = new u32[n_expanded];
    u32 *result2 = new u32[n_expanded];

    for (u32 i = 0; i < n_expanded; ++i)
    {
        result1[i] = montMod.to_u32(a_mont1[i]);
        result2[i] = montMod.to_u32(a_mont_simd_result[i]);
    }

    // 比较结果
    bool equal = compare_arrays(result1, result2, n_expanded);

    if (equal)
    {
        printf("测试通过: 常规NTT逆变换和SIMD NTT逆变换结果一致\n");
    }
    else
    {
        printf("测试失败: 常规NTT逆变换和SIMD NTT逆变换结果不一致\n");

        // 打印部分数据进行调试
        printf("前10个元素比较:\n");
        for (u32 i = 0; i < std::min(10u, n_expanded); ++i)
        {
            printf("a[%u]: 常规=%u, SIMD=%u\n", i, result1[i], result2[i]);
        }
    }

    // 释放内存
    delete[] a;
    delete[] a_mont1;
    delete[] a_mont2;
    delete[] a_mont2_copy;
    delete[] a_mont_simd;
    delete[] a_mont_simd_result;
    delete[] result1;
    delete[] result2;
}

// 测试NTT正变换后接逆变换的正确性
void test_ntt_forward_inverse_roundtrip(u32 n, u32 p, u32_mont omega_mont)
{
    printf("测试 n=%u, p=%u 的NTT正逆变换循环的正确性\n", n, p);

    // 确保n是2的幂
    u32 n_expanded = expand_n(n);

    // 生成随机数据
    u32 *original_data = new u32[n_expanded];
    generate_random_data(original_data, n_expanded, p);

    // 创建Montgomery域的数据
    MontMod montMod(p);
    u32_mont *a_mont = new u32_mont[n_expanded];

    // 转换为Montgomery域
    for (u32 i = 0; i < n_expanded; ++i)
    {
        a_mont[i] = montMod.from_u32(original_data[i]);
    }

    // 对输入数据执行bit-reverse排列
    u32_mont *a_mont_bit_reversed = new u32_mont[n_expanded];
    std::copy(a_mont, a_mont + n_expanded, a_mont_bit_reversed);
    bit_reverse_permute(a_mont_bit_reversed, n_expanded);

    // 转换为SIMD格式
    u32x4_mont *a_mont_simd = new u32x4_mont[n_expanded / 4];
    for (u32 i = 0; i < n_expanded / 4; ++i)
    {
        a_mont_simd[i] = vdupq_n_u32(0);
        for (u32 j = 0; j < 4; ++j)
        {
            a_mont_simd[i] = set_lane(a_mont_simd[i], a_mont_bit_reversed[i * 4 + j], j);
        }
    }

    // 执行正变换
    ntt_forward_mont_simd(a_mont_simd, n_expanded, p, omega_mont);

    // 执行逆变换
    ntt_inverse_dit_mont_simd(a_mont_simd, n_expanded, p, omega_mont);

    // 将结果转换回普通数组
    u32_mont *result_mont = new u32_mont[n_expanded];
    for (u32 i = 0; i < n_expanded / 4; ++i)
    {
        for (u32 j = 0; j < 4; ++j)
        {
            result_mont[i * 4 + j] = get_lane(a_mont_simd[i], j);
        }
    }

    // 对结果执行bit-reverse排列，恢复原始顺序
    bit_reverse_permute(result_mont, n_expanded);

    // 转换回普通域
    u32 *result = new u32[n_expanded];
    for (u32 i = 0; i < n_expanded; ++i)
    {
        result[i] = montMod.to_u32(result_mont[i]);
    }

    // 比较原始数据和结果
    bool equal = compare_arrays(original_data, result, n_expanded);

    if (equal)
    {
        printf("测试通过: 正变换后接逆变换恢复了原始数据\n");
    }
    else
    {
        printf("测试失败: 正变换后接逆变换未恢复原始数据\n");

        // 打印部分数据进行调试
        printf("前10个元素比较:\n");
        for (u32 i = 0; i < std::min(10u, n_expanded); ++i)
        {
            printf("a[%u]: 原始=%u, 结果=%u\n", i, original_data[i], result[i]);
        }
    }

    // 释放内存
    delete[] original_data;
    delete[] a_mont;
    delete[] a_mont_bit_reversed;
    delete[] a_mont_simd;
    delete[] result_mont;
    delete[] result;
}

int main()
{
    // 设置素数p和原根omega
    u32 p = 257; // 一个简单的素数，可以根据需要更改
    MontMod montMod(p);
    u32_mont omega_mont = montMod.from_u32(3); // 假设3是mod p的原根

    // 测试不同大小的输入
    test_ntt_inverse_consistency(8, p, omega_mont);
    test_ntt_inverse_consistency(16, p, omega_mont);
    test_ntt_inverse_consistency(32, p, omega_mont);
    test_ntt_inverse_consistency(64, p, omega_mont);

    printf("\n");

    // 测试正逆变换的循环正确性
    test_ntt_forward_inverse_roundtrip(8, p, omega_mont);
    test_ntt_forward_inverse_roundtrip(16, p, omega_mont);
    test_ntt_forward_inverse_roundtrip(32, p, omega_mont);
    test_ntt_forward_inverse_roundtrip(64, p, omega_mont);

    // 如果系统资源允许，可以测试更大的输入
    // test_ntt_inverse_consistency(128, p, omega_mont);
    // test_ntt_inverse_consistency(256, p, omega_mont);
    // test_ntt_forward_inverse_roundtrip(128, p, omega_mont);
    // test_ntt_forward_inverse_roundtrip(256, p, omega_mont);

    return 0;
}