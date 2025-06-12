#include "src/include/ntt.h"
#include "src/include/OpenMP/ntt.h"
#include "src/include/CRT/ntt.h"
#include "src/include/Barrett/ntt.h"
#include "src/include/OpenMP_Barrett/ntt.h"
#include "src/include/MPI/ntt.h"

#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <sys/time.h>
#include <unistd.h>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "src/include/MPI/config.h"

std::string path_read = "/home/hexay/projects/Lab Proj/Parallel-Programing-Project/.nttdata/";
std::string path_write = "/home/hexay/projects/Lab Proj/Parallel-Programing-Project/.results/";

template <typename T>
void fRead(T *a, T *b, T *n, T *p, T input_id)
{
    std::string str2 = std::to_string(input_id);
    std::string strin = path_read + str2 + ".in";
    std::ifstream fin(strin);
    fin >> *n >> *p;
    for (T i = 0; i < *n; i++)
    {
        fin >> a[i];
    }
    for (T i = 0; i < *n; i++)
    {
        fin >> b[i];
    }
}

// 线程开销测试
void test_thread_overhead()
{
    const int iterations = 1000;
    const int array_size = 131072;
    std::vector<double> data(array_size, 1.0);

    std::cout << "\n=== 线程创建开销测试 ===" << std::endl;

    // 测试不同线程数的开销
    for (int num_threads = 1; num_threads <= 8; num_threads *= 2)
    {
#ifdef _OPENMP
        omp_set_num_threads(num_threads);
#endif

        auto start = std::chrono::high_resolution_clock::now();

        for (int iter = 0; iter < iterations; iter++)
        {
            double sum = 0.0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : sum)
#endif
            for (int i = 0; i < array_size; i++)
            {
                sum += data[i] * i;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        std::cout << "线程数: " << num_threads
                  << ", 平均时间: " << duration.count() / (double)iterations << " μs" << std::endl;
    }
}

// 内存带宽竞争测试
void test_memory_bandwidth()
{
    const int array_size = 131072;
    const int iterations = 100;
    std::vector<u64> src(array_size);
    std::vector<u64> dst(array_size);

    for (int i = 0; i < array_size; i++)
    {
        src[i] = i;
    }

    std::cout << "\n=== 内存带宽竞争测试 ===" << std::endl;

    for (int num_threads = 1; num_threads <= 8; num_threads *= 2)
    {
#ifdef _OPENMP
        omp_set_num_threads(num_threads);
#endif

        auto start = std::chrono::high_resolution_clock::now();

        for (int iter = 0; iter < iterations; iter++)
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < array_size; i++)
            {
                dst[i] = src[i] * 2 + 1;
            }
        }

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        double bandwidth = (double)(array_size * sizeof(u64) * 2 * iterations) /
                           (duration.count() / 1e6) / (1024 * 1024 * 1024);

        std::cout << "线程数: " << num_threads
                  << ", 平均时间: " << duration.count() / (double)iterations << " μs"
                  << ", 带宽: " << bandwidth << " GB/s" << std::endl;
    }
}

// NTT性能详细分析
void detailed_ntt_performance_analysis()
{
    u64 a[300000], b[300000], ab[300000];
    u64 n_, p_;

    // 使用大规模数据（input_id = 2）
    fRead<u64>(a, b, &n_, &p_, 2);

    std::cout << "\n=== NTT性能详细分析 (n=" << n_ << ", p=" << p_ << ") ===" << std::endl;

    // 测试不同OpenMP线程数的性能
    for (int num_threads = 1; num_threads <= 8; num_threads *= 2)
    {
#ifdef _OPENMP
        omp_set_num_threads(num_threads);
#endif

        const int test_runs = 5;
        double total_time = 0.0;

        for (int run = 0; run < test_runs; run++)
        {
            memset(ab, 0, sizeof(ab));

            auto start = std::chrono::high_resolution_clock::now();
            poly_multiply_ntt_omp_Barrett(a, b, ab, n_, p_);
            auto end = std::chrono::high_resolution_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_time += duration.count();
        }

        double avg_time = total_time / test_runs;
        double speedup = (num_threads == 1) ? 1.0 : (avg_time / avg_time); // Will be updated
        static double single_thread_time = 0.0;

        if (num_threads == 1)
        {
            single_thread_time = avg_time;
            speedup = 1.0;
        }
        else
        {
            speedup = single_thread_time / avg_time;
        }

        double efficiency = speedup / num_threads * 100.0;

        std::cout << "线程数: " << num_threads
                  << ", 平均时间: " << std::fixed << std::setprecision(2) << avg_time << " μs"
                  << ", 加速比: " << std::setprecision(2) << speedup
                  << ", 效率: " << std::setprecision(1) << efficiency << "%" << std::endl;
    }
}

// 任务粒度分析
void analyze_task_granularity()
{
    std::cout << "\n=== 任务粒度分析 ===" << std::endl;

    u64 n = 131072;
    std::cout << "问题规模 n = " << n << std::endl;
    std::cout << "NTT总轮数: " << (int)log2(n) << std::endl;

    std::cout << "\n各轮次可并行任务数量:" << std::endl;
    for (u64 mid = 1; mid < n; mid <<= 1)
    {
        u64 num_tasks = n / (mid << 1);
        u64 work_per_task = mid;
        std::cout << "轮次 " << (int)log2(mid) + 1
                  << ": mid=" << mid
                  << ", 任务数=" << num_tasks
                  << ", 每任务工作量=" << work_per_task << std::endl;
    }
}

// 性能热力图数据收集
void collect_performance_heatmap_data()
{
    std::cout << "\n=== 收集性能热力图数据 ===" << std::endl;

    u64 a[300000], b[300000], ab[300000];
    u64 n_, p_;

    // 使用大规模数据
    fRead<u64>(a, b, &n_, &p_, 2);

    std::ofstream outfile(path_write + "performance_heatmap.csv");
    outfile << "threads,avg_time_us,speedup,efficiency\n";

    double single_thread_time = 0.0;

    for (int num_threads = 1; num_threads <= 8; num_threads++)
    {
#ifdef _OPENMP
        omp_set_num_threads(num_threads);
#endif

        const int test_runs = 3;
        double total_time = 0.0;

        for (int run = 0; run < test_runs; run++)
        {
            memset(ab, 0, sizeof(ab));

            auto start = std::chrono::high_resolution_clock::now();
            poly_multiply_ntt_omp_Barrett(a, b, ab, n_, p_);
            auto end = std::chrono::high_resolution_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            total_time += duration.count();
        }

        double avg_time = total_time / test_runs;

        if (num_threads == 1)
        {
            single_thread_time = avg_time;
        }

        double speedup = single_thread_time / avg_time;
        double efficiency = speedup / num_threads * 100.0;

        outfile << num_threads << "," << avg_time << "," << speedup << "," << efficiency << "\n";

        std::cout << "线程数 " << num_threads << " 完成" << std::endl;
    }

    outfile.close();
}

int main(int argc, char *argv[])
{
    MPI_INIT(argc, argv);

    MPI_GET_RANK(rank);

    MPI_ONLY_MAIN
    {
        std::cout << "=== 混合并行性能退化分析实验 ===" << std::endl;

        // 确保输出目录存在
        system("mkdir -p /home/hexay/projects/Lab\\ Proj/Parallel-Programing-Project/.results/");

        // 运行各项测试
        analyze_task_granularity();
        test_thread_overhead();
        test_memory_bandwidth();
        detailed_ntt_performance_analysis();
        collect_performance_heatmap_data();

        std::cout << "\n=== 实验完成 ===" << std::endl;
    }

    MPI_FINALIZE();
    return 0;
}