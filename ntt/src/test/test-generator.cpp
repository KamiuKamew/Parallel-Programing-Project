#include "../include/ntt.h"
#include "../include/simd/ntt.h"

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
              << "  -simd           执行SIMD优化的NTT乘法 (默认)\n"
              << "  -all            执行所有算法并比较结果\n"
              << "  -verify         验证结果的正确性\n"
              << "  -h, --help      显示此帮助信息\n";
}

int main(int argc, char *argv[])
{
    int n = 32768;          // 默认多项式大小
    int iterations = 1;     // 默认迭代次数
    int p = 998244353;      // 默认模数
    bool run_naive = false; // 是否运行朴素算法
    bool run_ntt = false;   // 是否运行普通NTT
    bool run_simd = true;   // 默认运行SIMD NTT
    bool verify = false;    // 是否验证结果

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
        else if (strcmp(argv[i], "-simd") == 0)
        {
            run_simd = true;
        }
        else if (strcmp(argv[i], "-all") == 0)
        {
            run_naive = run_ntt = run_simd = true;
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
    int *a = new int[n_expanded];
    int *b = new int[n_expanded];
    int *ab_naive = new int[n_expanded];
    int *ab_ntt = new int[n_expanded];
    int *ab_ntt_simd = new int[n_expanded];

    // 初始化结果数组
    std::memset(ab_naive, 0, n_expanded * sizeof(int));
    std::memset(ab_ntt, 0, n_expanded * sizeof(int));
    std::memset(ab_ntt_simd, 0, n_expanded * sizeof(int));

    // 生成随机多项式
    generate_random_polynomial(a, n, p);
    generate_random_polynomial(b, n, p);

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
            poly_multiply_naive(a, b, ab_naive, n, p);
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        std::cout << "朴素乘法: " << duration / 1000.0 << " ms (平均每次 "
                  << duration / (double)iterations / 1000.0 << " ms)\n";
    }

    // 测试NTT乘法
    if (run_ntt)
    {
        auto start = std::chrono::high_resolution_clock::now();

        for (int iter = 0; iter < iterations; ++iter)
        {
            std::memset(ab_ntt, 0, n_expanded * sizeof(int));
            poly_multiply_ntt(a, b, ab_ntt, n, p);
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        std::cout << "NTT乘法: " << duration / 1000.0 << " ms (平均每次 "
                  << duration / (double)iterations / 1000.0 << " ms)\n";

        // 验证与朴素乘法结果
        if (verify && run_naive)
        {
            bool valid = verify_results(ab_naive, ab_ntt, n_expanded);
            std::cout << "NTT结果验证: " << (valid ? "正确" : "不正确") << "\n";
        }
    }

    // 测试SIMD NTT乘法
    if (run_simd)
    {
        auto start = std::chrono::high_resolution_clock::now();

        for (int iter = 0; iter < iterations; ++iter)
        {
            std::memset(ab_ntt_simd, 0, n_expanded * sizeof(int));
            poly_multiply_ntt_simd(a, b, ab_ntt_simd, n, p);
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        std::cout << "SIMD NTT乘法: " << duration / 1000.0 << " ms (平均每次 "
                  << duration / (double)iterations / 1000.0 << " ms)\n";

        // 验证与朴素乘法结果
        if (verify && run_naive)
        {
            bool valid = verify_results(ab_naive, ab_ntt_simd, n_expanded);
            std::cout << "SIMD NTT结果验证: " << (valid ? "正确" : "不正确") << "\n";
        }

        // 验证与普通NTT结果
        if (verify && run_ntt)
        {
            bool valid = verify_results(ab_ntt, ab_ntt_simd, n_expanded);
            std::cout << "SIMD vs NTT结果验证: " << (valid ? "正确" : "不正确") << "\n";
        }
    }

    // 如果使用了所有算法，打印结果以便比较
    if (run_naive && run_ntt && run_simd && n <= 16)
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

        std::cout << "SIMD: ";
        for (int i = 0; i < std::min(n_expanded, 16); ++i)
        {
            std::cout << ab_ntt_simd[i] << " ";
        }
        std::cout << "\n";
    }

    // 释放内存
    delete[] a;
    delete[] b;
    delete[] ab_naive;
    delete[] ab_ntt;
    delete[] ab_ntt_simd;

    return 0;
}