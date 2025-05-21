#include "../include/ntt.h"
#include "../include/CRT/ntt.h"

#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cstring>

// 生成随机多项式
void generate_random_polynomial(int *poly, int n, int p)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> dis(0, p - 1);

    for (int i = 0; i < n; ++i)
    {
        poly[i] = dis(gen);
    }

    // 填充剩余部分为0
    int n_expanded = 1;
    while (n_expanded < 2 * n - 1)
        n_expanded <<= 1;

    for (int i = n; i < n_expanded; ++i)
    {
        poly[i] = 0;
    }
}

// 验证结果是否一致
bool verify_results(int *result1, int *result2, int n)
{
    for (int i = 0; i < n; ++i)
    {
        if (result1[i] != result2[i])
        {
            return false;
        }
    }
    return true;
}

// 打印使用说明
void print_usage(const char *program_name)
{
    std::cerr << "使用方法: " << program_name << " [选项]\n"
              << "选项:\n"
              << "  -n <size>       多项式的大小 (默认: 4)\n"
              << "  -iter <num>     测试迭代次数 (默认: 100)\n"
              << "  -p <prime>      使用的模数 (默认: 998244353)\n"
              << "  -naive          执行朴素乘法\n"
              << "  -ntt            执行NTT乘法\n"
              << "  -crt           执行crt优化的NTT乘法 (默认)\n"
              << "  -all            执行所有算法并比较结果\n"
              << "  -verify         验证结果的正确性\n"
              << "  -h, --help      显示此帮助信息\n";
}

int main(int argc, char *argv[])
{
    int n = 1024;           // 默认多项式大小
    int iterations = 1;     // 默认迭代次数
    int p = 7340033;        // 默认模数
    bool run_naive = false; // 是否运行朴素算法
    bool run_ntt = true;    // 是否运行普通NTT
    bool run_crt = false;   // 默认运行crt NTT
    bool verify = true;     // 是否验证结果

    // 解析命令行参数
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "-n") == 0 && i + 1 < argc)
        {
            n = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-iter") == 0 && i + 1 < argc)
        {
            iterations = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-p") == 0 && i + 1 < argc)
        {
            p = std::atoi(argv[++i]);
        }
        else if (strcmp(argv[i], "-naive") == 0)
        {
            run_naive = true;
        }
        else if (strcmp(argv[i], "-ntt") == 0)
        {
            run_ntt = true;
        }
        else if (strcmp(argv[i], "-crt") == 0)
        {
            run_crt = true;
        }
        else if (strcmp(argv[i], "-all") == 0)
        {
            run_naive = run_ntt = run_crt = true;
        }
        else if (strcmp(argv[i], "-verify") == 0)
        {
            verify = true;
        }
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
        {
            print_usage(argv[0]);
            return 0;
        }
        else
        {
            std::cerr << "未知选项: " << argv[i] << std::endl;
            print_usage(argv[0]);
            return 1;
        }
    }

    // 计算需要的扩展数组大小
    int n_expanded = 1;
    while (n_expanded < 2 * n - 1)
        n_expanded <<= 1;

    // 分配内存
    int *a_int = new int[n_expanded];
    int *b_int = new int[n_expanded];
    int *ab_naive = new int[n_expanded];
    int *ab_ntt = new int[n_expanded];
    u64 *a_crt = new u64[n_expanded];
    u64 *b_crt = new u64[n_expanded];
    u64 *ab_ntt_crt_u64 = new u64[n_expanded];
    int *ab_ntt_crt_int = new int[n_expanded];

    // 初始化结果数组
    std::memset(ab_naive, 0, n_expanded * sizeof(int));
    std::memset(ab_ntt, 0, n_expanded * sizeof(int));
    std::memset(ab_ntt_crt_u64, 0, n_expanded * sizeof(u64));
    std::memset(ab_ntt_crt_int, 0, n_expanded * sizeof(int));

    // 生成随机多项式
    generate_random_polynomial(a_int, n, p);
    generate_random_polynomial(b_int, n, p);

    // Populate u64 versions for CRT
    for (int i = 0; i < n_expanded; ++i)
    {
        a_crt[i] = static_cast<u64>(a_int[i]);
        b_crt[i] = static_cast<u64>(b_int[i]);
    }

    std::cout << "测试配置:\n";
    std::cout << "- 多项式大小: " << n << "\n";
    std::cout << "- 迭代次数: " << iterations << "\n";
    std::cout << "- 模数: " << p << "\n";

    // 测试朴素乘法
    if (run_naive)
    {
        auto start = std::chrono::high_resolution_clock::now();

        for (int iter = 0; iter < iterations; ++iter)
        {
            std::memset(ab_naive, 0, n_expanded * sizeof(int));
            poly_multiply_naive(a_int, b_int, ab_naive, n, p);
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        std::cout << "朴素乘法: " << duration / 1000.0 << " ms (平均每次 "
                  << duration / (double)iterations / 1000.0 << " ms)\n";

        std::cout << "朴素乘法结果: ";
        for (int i = 0; i < n_expanded; ++i)
        {
            std::cout << ab_naive[i] << " ";
        }
        std::cout << std::endl;
    }

    // 测试NTT乘法
    if (run_ntt)
    {
        auto start = std::chrono::high_resolution_clock::now();

        for (int iter = 0; iter < iterations; ++iter)
        {
            std::memset(ab_ntt, 0, n_expanded * sizeof(int));
            poly_multiply_ntt(a_int, b_int, ab_ntt, n, p);
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        std::cout << "NTT乘法: " << duration / 1000.0 << " ms (平均每次 "
                  << duration / (double)iterations / 1000.0 << " ms)\n";

        std::cout << "NTT乘法结果: ";
        for (int i = 0; i < n_expanded; ++i)
        {
            std::cout << ab_ntt[i] << " ";
        }
        std::cout << std::endl;

        // 验证与朴素乘法结果
        if (verify && run_naive)
        {
            bool valid = verify_results(ab_naive, ab_ntt, n_expanded);
            std::cout << "NTT结果验证: " << (valid ? "正确" : "不正确") << "\n";
        }
    }

    // 测试crt NTT乘法
    if (run_crt)
    {
        auto start = std::chrono::high_resolution_clock::now();

        for (int iter = 0; iter < iterations; ++iter)
        {
            std::memset(ab_ntt_crt_u64, 0, n_expanded * sizeof(u64));
            poly_multiply_ntt_crt(a_crt, b_crt, ab_ntt_crt_u64, n, p);
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        // Convert u64 result to int for verification and printing
        for (int i = 0; i < n_expanded; ++i)
        {
            ab_ntt_crt_int[i] = static_cast<int>(ab_ntt_crt_u64[i]);
        }

        std::cout << "crt NTT乘法: " << duration / 1000.0 << " ms (平均每次 "
                  << duration / (double)iterations / 1000.0 << " ms)\n";

        std::cout << "crt NTT乘法结果: ";
        for (int i = 0; i < n_expanded; ++i)
        {
            std::cout << ab_ntt_crt_int[i] << " ";
        }
        std::cout << std::endl;

        // 验证与朴素乘法结果
        if (verify && run_naive)
        {
            bool valid = verify_results(ab_naive, ab_ntt_crt_int, n_expanded);
            std::cout << "crt NTT结果验证: " << (valid ? "正确" : "不正确") << "\n";
        }

        // 验证与普通NTT结果
        if (verify && run_ntt)
        {
            bool valid = verify_results(ab_ntt, ab_ntt_crt_int, n_expanded);
            std::cout << "crt vs NTT结果验证: " << (valid ? "正确" : "不正确") << "\n";
        }
    }

    // 如果使用了所有算法，打印结果以便比较
    if (run_naive && run_ntt && run_crt && n <= 16)
    {
        std::cout << "\n结果比较 (前" << std::min(n_expanded, 16) << "个元素):\n";

        std::cout << "朴素: ";
        for (int i = 0; i < std::min(n_expanded, 16); ++i)
        {
            std::cout << ab_naive[i] << " ";
        }
        std::cout << "\n";

        std::cout << "NTT:  ";
        for (int i = 0; i < std::min(n_expanded, 16); ++i)
        {
            std::cout << ab_ntt[i] << " ";
        }
        std::cout << "\n";

        std::cout << "crt: ";
        for (int i = 0; i < std::min(n_expanded, 16); ++i)
        {
            std::cout << ab_ntt_crt_int[i] << " ";
        }
        std::cout << "\n";
    }

    // 释放内存
    delete[] a_int;
    delete[] b_int;
    delete[] ab_naive;
    delete[] ab_ntt;
    delete[] a_crt;
    delete[] b_crt;
    delete[] ab_ntt_crt_u64;
    delete[] ab_ntt_crt_int;

    return 0;
}