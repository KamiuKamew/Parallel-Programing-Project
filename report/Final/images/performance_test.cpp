#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <iomanip>
#include <fstream>
#include <string>
#include <algorithm>
#include <thread>
#include <mutex>
#include <map>
#include <mpi.h>
#include <omp.h>

#include "src/include/unified/unified.h"
#include "src/include/general/utils.h"

// 测试配置
struct TestConfig
{
    std::vector<u32> sizes = {1024, 8192, 131072}; // 减少到3个典型规模
    u32 p = 998244353;                             // 标准模数
    u32 omega = 3;                                 // 原根
    int iterations = 5;                            // 减少到5次迭代
    int warmup_runs = 2;                           // 减少预热次数
};

// 测试结果结构
struct TestResult
{
    ModMethod mod_method;
    ParallelMethod parallel_method;
    u32 size;
    double avg_time_ms;
    double min_time_ms;
    double max_time_ms;
    double std_dev_ms;
    bool correctness;
};

// 全局结果存储
std::vector<TestResult> all_results;
std::mutex results_mutex;

// 生成随机测试数据
void generate_test_data(u32 *a, u32 *b, u32 n, u32 p)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<u32> dis(0, p - 1);

    for (u32 i = 0; i < n; ++i)
    {
        a[i] = dis(gen);
        b[i] = dis(gen);
    }
}

// 验证结果正确性
bool verify_result(u32 *expected, u32 *actual, u32 n, u32 p)
{
    (void)p; // 消除未使用参数警告
    for (u32 i = 0; i < 2 * n - 1; ++i)
    {
        if (expected[i] != actual[i])
        {
            std::cerr << "Correctness check failed at index " << i
                      << ": expected " << expected[i]
                      << ", got " << actual[i] << std::endl;
            return false;
        }
    }
    return true;
}

// 运行单个测试
TestResult run_single_test(ModMethod mod_method, ParallelMethod parallel_method,
                           u32 size, const TestConfig &config)
{
    u32 n = size;
    u32 p = config.p;
    u32 omega = config.omega;

    // 分配内存
    u32 *a = new u32[n];
    u32 *b = new u32[n];
    u32 *result = new u32[2 * n - 1];
    u32 *expected = new u32[2 * n - 1];

    // 生成测试数据
    generate_test_data(a, b, n, p);

    // 计算期望结果（使用串行Montgomery作为基准）
    {
        u32 *a_copy = new u32[n];
        u32 *b_copy = new u32[n];
        std::copy(a, a + n, a_copy);
        std::copy(b, b + n, b_copy);

        poly_multiply_unified_choosewith(a_copy, b_copy, expected, n, p, omega,
                                         ModMethod::MONTGOMERY, ParallelMethod::SERIAL);

        delete[] a_copy;
        delete[] b_copy;
    }

    // 预热运行
    for (int warmup = 0; warmup < config.warmup_runs; ++warmup)
    {
        std::copy(a, a + n, result);
        poly_multiply_unified_choosewith(a, b, result, n, p, omega,
                                         mod_method, parallel_method);
    }

    // 性能测试
    std::vector<double> times;
    bool correctness = true;

    for (int iter = 0; iter < config.iterations; ++iter)
    {
        // 重置数据
        std::copy(a, a + n, result);

        auto start = std::chrono::high_resolution_clock::now();

        poly_multiply_unified_choosewith(a, b, result, n, p, omega,
                                         mod_method, parallel_method);

        auto end = std::chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        double time_ms = duration.count() / 1000.0;
        times.push_back(time_ms);

        // 验证正确性（只在第一次迭代时验证）
        if (iter == 0)
        {
            correctness = verify_result(expected, result, n, p);
        }
    }

    // 计算统计信息
    double sum = 0.0;
    for (double time : times)
    {
        sum += time;
    }
    double avg_time = sum / times.size();

    double min_time = *std::min_element(times.begin(), times.end());
    double max_time = *std::max_element(times.begin(), times.end());

    // 计算标准差
    double variance = 0.0;
    for (double time : times)
    {
        variance += (time - avg_time) * (time - avg_time);
    }
    double std_dev = std::sqrt(variance / times.size());

    // 清理内存
    delete[] a;
    delete[] b;
    delete[] result;
    delete[] expected;

    return {mod_method, parallel_method, size, avg_time, min_time, max_time, std_dev, correctness};
}

// 运行所有测试
void run_all_tests(const TestConfig &config, int rank = 0)
{
    std::vector<ModMethod> mod_methods = {
        ModMethod::NAIVE,
        ModMethod::MONTGOMERY,
        ModMethod::BARRETT};

    std::vector<ParallelMethod> parallel_methods = {
        ParallelMethod::SERIAL,
        ParallelMethod::OPENMP,
        ParallelMethod::MPI,
        ParallelMethod::GPU};

    int total_tests = mod_methods.size() * parallel_methods.size() * config.sizes.size();
    int current_test = 0;

    if (rank == 0)
    {
        std::cout << "开始性能测试..." << std::endl;
        std::cout << "总测试数: " << total_tests << std::endl;
        std::cout << "每个测试重复次数: " << config.iterations << std::endl;
        std::cout << std::endl;
    }

    for (auto mod_method : mod_methods)
    {
        for (auto parallel_method : parallel_methods)
        {
            for (auto size : config.sizes)
            {
                current_test++;

                // 只有主进程输出进度信息
                if (rank == 0)
                {
                    std::cout << "测试 " << current_test << "/" << total_tests
                              << " - 算法: " << static_cast<int>(mod_method)
                              << ", 并行策略: " << static_cast<int>(parallel_method)
                              << ", 规模: " << size << std::endl;
                }

                TestResult result = run_single_test(mod_method, parallel_method, size, config);

                {
                    std::lock_guard<std::mutex> lock(results_mutex);
                    all_results.push_back(result);
                }

                // 只有主进程输出结果
                if (rank == 0)
                {
                    std::cout << "  平均时间: " << std::fixed << std::setprecision(3)
                              << result.avg_time_ms << " ms" << std::endl;
                    std::cout << "  正确性: " << (result.correctness ? "通过" : "失败") << std::endl;
                    std::cout << std::endl;
                }
            }
        }
    }
}

// 输出结果到CSV文件
void output_results_to_csv(const std::string &filename)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "无法创建文件: " << filename << std::endl;
        return;
    }

    // CSV头部
    file << "Algorithm,Parallel_Strategy,Size,Avg_Time_ms,Min_Time_ms,Max_Time_ms,Std_Dev_ms,Correctness\n";

    // 算法名称映射
    std::vector<std::string> algorithm_names = {"Naive", "Montgomery", "Barrett"};
    std::vector<std::string> parallel_names = {"Serial", "OpenMP", "MPI", "GPU"};

    // 创建并行方法到索引的映射
    std::map<ParallelMethod, int> parallel_to_index = {
        {ParallelMethod::SERIAL, 0},
        {ParallelMethod::OPENMP, 1},
        {ParallelMethod::MPI, 2},
        {ParallelMethod::GPU, 3}};

    for (const auto &result : all_results)
    {
        int parallel_index = parallel_to_index[result.parallel_method];
        file << algorithm_names[static_cast<int>(result.mod_method)] << ","
             << parallel_names[parallel_index] << ","
             << result.size << ","
             << std::fixed << std::setprecision(3) << result.avg_time_ms << ","
             << result.min_time_ms << ","
             << result.max_time_ms << ","
             << result.std_dev_ms << ","
             << (result.correctness ? "Pass" : "Fail") << "\n";
    }

    file.close();
    std::cout << "结果已保存到: " << filename << std::endl;
}

// 输出性能总结
void output_performance_summary()
{
    std::cout << "\n=== 性能测试总结 ===" << std::endl;

    std::vector<std::string> algorithm_names = {"Naive", "Montgomery", "Barrett"};
    std::vector<std::string> parallel_names = {"Serial", "OpenMP", "MPI", "GPU"};

    // 创建并行方法到索引的映射
    std::map<ParallelMethod, int> parallel_to_index = {
        {ParallelMethod::SERIAL, 0},
        {ParallelMethod::OPENMP, 1},
        {ParallelMethod::MPI, 2},
        {ParallelMethod::GPU, 3}};

    for (size_t i = 0; i < algorithm_names.size(); ++i)
    {
        for (size_t j = 0; j < parallel_names.size(); ++j)
        {
            std::cout << "\n"
                      << algorithm_names[i] << " + " << parallel_names[j] << ":" << std::endl;

            for (const auto &result : all_results)
            {
                if (static_cast<int>(result.mod_method) == static_cast<int>(i) &&
                    parallel_to_index[result.parallel_method] == static_cast<int>(j))
                {
                    std::cout << "  规模 " << result.size << ": "
                              << std::fixed << std::setprecision(3) << result.avg_time_ms
                              << " ms (" << (result.correctness ? "正确" : "错误") << ")" << std::endl;
                }
            }
        }
    }
}

int main()
{
    // 初始化MPI
    int rank, size;
    MPI_Init(nullptr, nullptr);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 设置OpenMP线程数
    omp_set_num_threads(8); // 使用8个线程，避免过度并行化

    // 只有主进程输出信息
    if (rank == 0)
    {
        TestConfig config;

        std::cout << "NTT统一框架性能测试" << std::endl;
        std::cout << "==================" << std::endl;
        std::cout << "测试规模: ";
        for (size_t i = 0; i < config.sizes.size(); ++i)
        {
            std::cout << config.sizes[i];
            if (i < config.sizes.size() - 1)
                std::cout << ", ";
        }
        std::cout << std::endl;
        std::cout << "模数: " << config.p << std::endl;
        std::cout << "原根: " << config.omega << std::endl;
        std::cout << "每个组合重复次数: " << config.iterations << std::endl;
        std::cout << std::endl;

        // 运行测试
        run_all_tests(config, rank);

        // 输出结果
        output_results_to_csv("performance_results.csv");
        output_performance_summary();
    }
    else
    {
        // 非主进程等待
        TestConfig config;
        run_all_tests(config, rank);
    }

    // 清理MPI
    MPI_Finalize();

    return 0;
}