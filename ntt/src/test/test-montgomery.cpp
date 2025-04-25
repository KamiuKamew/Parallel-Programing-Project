#include "../include/op.h"

#include <iostream>
#include <cassert>
#include <random>

// 打印测试结果
void print_success(const std::string &test_name)
{
    std::cout << "[✓] " << test_name << " 测试通过" << std::endl;
}

// 检查Montgomery模运算和普通模运算的结果是否一致
template <typename T>
void assert_equal(T expected, T actual, const std::string &operation)
{
    if (expected != actual)
    {
        std::cout << "[✗] " << operation << " 测试失败! 预期: " << expected << ", 实际: " << actual << std::endl;
        assert(false);
    }
}

int main()
{
    // 使用常用的NTT模数进行测试
    u32 mod = 998244353; // 2^23 * 119 + 1

    // 创建普通模运算和Montgomery模运算对象
    Mod normal_mod(mod);
    Montgomery montgomery(mod);

    // 测试基本运算
    std::cout << "===== 测试基本运算 =====" << std::endl;

    // 测试 from_u32 和 to_u32
    {
        u32 a = 123456789;
        u32 mont_a = montgomery.from_u32(a);
        u32 back_a = montgomery.to_u32(mont_a);
        assert_equal(a, back_a, "转换测试 (a -> Montgomery域 -> a)");
        print_success("转换测试");
    }

    // 测试加法
    {
        u32 a = 123456789;
        u32 b = 987654321;

        u32 normal_add = normal_mod.add(a, b);

        u32 mont_a = montgomery.from_u32(a);
        u32 mont_b = montgomery.from_u32(b);
        u32 mont_add = montgomery.add(mont_a, mont_b);
        u32 result_add = montgomery.to_u32(mont_add);

        assert_equal(normal_add, result_add, "加法测试");
        print_success("加法测试");
    }

    // 测试减法
    {
        u32 a = 987654321;
        u32 b = 123456789;

        u32 normal_sub = normal_mod.sub(a, b);

        u32 mont_a = montgomery.from_u32(a);
        u32 mont_b = montgomery.from_u32(b);
        u32 mont_sub = montgomery.sub(mont_a, mont_b);
        u32 result_sub = montgomery.to_u32(mont_sub);

        assert_equal(normal_sub, result_sub, "减法测试");
        print_success("减法测试");
    }

    // 测试乘法
    {
        u32 a = 123456789;
        u32 b = 987654321;

        u32 normal_mul = normal_mod.mul(a, b);

        u32 mont_a = montgomery.from_u32(a);
        u32 mont_b = montgomery.from_u32(b);
        u32 mont_mul = montgomery.mul(mont_a, mont_b);
        u32 result_mul = montgomery.to_u32(mont_mul);

        assert_equal(normal_mul, result_mul, "乘法测试");
        print_success("乘法测试");
    }

    // 测试幂运算
    {
        u32 base = 123456;
        u32 exp = 12345;

        u32 normal_pow = normal_mod.pow(base, exp);
        u32 mont_pow = montgomery.pow(base, exp);

        assert_equal(normal_pow, mont_pow, "幂运算测试");
        print_success("幂运算测试");
    }

    // 测试逆元计算
    {
        u32 a = 123456789;

        u32 normal_inv = normal_mod.inv(a);
        u32 mont_inv = montgomery.inv(a);

        assert_equal(normal_inv, mont_inv, "逆元计算测试");
        print_success("逆元计算测试");
    }

    // 测试复合运算
    std::cout << "\n===== 测试复合运算 =====" << std::endl;

    // 测试 (a * b + c) % mod
    {
        u32 a = 123456789;
        u32 b = 987654321;
        u32 c = 567890123;

        // 普通模运算计算 (a * b + c) % mod
        u32 normal_mul = normal_mod.mul(a, b);
        u32 normal_result = normal_mod.add(normal_mul, c);

        // Montgomery模运算计算 (a * b + c) % mod
        u32 mont_a = montgomery.from_u32(a);
        u32 mont_b = montgomery.from_u32(b);
        u32 mont_c = montgomery.from_u32(c);
        u32 mont_mul = montgomery.mul(mont_a, mont_b);
        u32 mont_result = montgomery.add(mont_mul, mont_c);
        u32 result = montgomery.to_u32(mont_result);

        assert_equal(normal_result, result, "复合运算测试 (a * b + c) % mod");
        print_success("复合运算 (a * b + c) % mod");
    }

    // 测试 (a * b - c) % mod
    {
        u32 a = 123456789;
        u32 b = 987654321;
        u32 c = 567890123;

        // 普通模运算计算 (a * b - c) % mod
        u32 normal_mul = normal_mod.mul(a, b);
        u32 normal_result = normal_mod.sub(normal_mul, c);

        // Montgomery模运算计算 (a * b - c) % mod
        u32 mont_a = montgomery.from_u32(a);
        u32 mont_b = montgomery.from_u32(b);
        u32 mont_c = montgomery.from_u32(c);
        u32 mont_mul = montgomery.mul(mont_a, mont_b);
        u32 mont_result = montgomery.sub(mont_mul, mont_c);
        u32 result = montgomery.to_u32(mont_result);

        assert_equal(normal_result, result, "复合运算测试 (a * b - c) % mod");
        print_success("复合运算 (a * b - c) % mod");
    }

    // 测试 (a * b * c) % mod
    {
        u32 a = 123456789;
        u32 b = 987654321;
        u32 c = 567890123;

        // 普通模运算计算 (a * b * c) % mod
        u32 normal_mul1 = normal_mod.mul(a, b);
        u32 normal_result = normal_mod.mul(normal_mul1, c);

        // Montgomery模运算计算 (a * b * c) % mod
        u32 mont_a = montgomery.from_u32(a);
        u32 mont_b = montgomery.from_u32(b);
        u32 mont_c = montgomery.from_u32(c);
        u32 mont_mul1 = montgomery.mul(mont_a, mont_b);
        u32 mont_result = montgomery.mul(mont_mul1, mont_c);
        u32 result = montgomery.to_u32(mont_result);

        assert_equal(normal_result, result, "复合运算测试 (a * b * c) % mod");
        print_success("复合运算 (a * b * c) % mod");
    }

    // 测试随机数据
    std::cout << "\n===== 测试随机数据 =====" << std::endl;

    // 创建随机数生成器
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<u32> dist(1, mod - 1);

    const int NUM_TESTS = 100;

    // 测试随机加法
    {
        int success_count = 0;
        for (int i = 0; i < NUM_TESTS; ++i)
        {
            u32 a = dist(gen);
            u32 b = dist(gen);

            u32 normal_add = normal_mod.add(a, b);

            u32 mont_a = montgomery.from_u32(a);
            u32 mont_b = montgomery.from_u32(b);
            u32 mont_add = montgomery.add(mont_a, mont_b);
            u32 result_add = montgomery.to_u32(mont_add);

            if (normal_add == result_add)
            {
                success_count++;
            }
        }
        std::cout << "随机加法测试: " << success_count << "/" << NUM_TESTS << " 通过" << std::endl;
        assert(success_count == NUM_TESTS);
    }

    // 测试随机乘法
    {
        int success_count = 0;
        for (int i = 0; i < NUM_TESTS; ++i)
        {
            u32 a = dist(gen);
            u32 b = dist(gen);

            u32 normal_mul = normal_mod.mul(a, b);

            u32 mont_a = montgomery.from_u32(a);
            u32 mont_b = montgomery.from_u32(b);
            u32 mont_mul = montgomery.mul(mont_a, mont_b);
            u32 result_mul = montgomery.to_u32(mont_mul);

            if (normal_mul == result_mul)
            {
                success_count++;
            }
        }
        std::cout << "随机乘法测试: " << success_count << "/" << NUM_TESTS << " 通过" << std::endl;
        assert(success_count == NUM_TESTS);
    }

    // 测试随机复合运算 (a * b + c * d) % mod
    {
        int success_count = 0;
        for (int i = 0; i < NUM_TESTS; ++i)
        {
            u32 a = dist(gen);
            u32 b = dist(gen);
            u32 c = dist(gen);
            u32 d = dist(gen);

            // 普通模运算
            u32 normal_mul1 = normal_mod.mul(a, b);
            u32 normal_mul2 = normal_mod.mul(c, d);
            u32 normal_result = normal_mod.add(normal_mul1, normal_mul2);

            // Montgomery模运算
            u32 mont_a = montgomery.from_u32(a);
            u32 mont_b = montgomery.from_u32(b);
            u32 mont_c = montgomery.from_u32(c);
            u32 mont_d = montgomery.from_u32(d);

            u32 mont_mul1 = montgomery.mul(mont_a, mont_b);
            u32 mont_mul2 = montgomery.mul(mont_c, mont_d);
            u32 mont_result = montgomery.add(mont_mul1, mont_mul2);
            u32 result = montgomery.to_u32(mont_result);

            if (normal_result == result)
            {
                success_count++;
            }
        }
        std::cout << "随机复合运算测试: " << success_count << "/" << NUM_TESTS << " 通过" << std::endl;
        assert(success_count == NUM_TESTS);
    }

    std::cout << "\n所有测试通过！Montgomery规约模运算实现正确。" << std::endl;

    return 0;
}
