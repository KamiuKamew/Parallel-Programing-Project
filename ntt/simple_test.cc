#include "src/include/OpenMP_Barrett/ntt.h"
#include <chrono>
#include <fstream>
#include <iostream>
#include <cstring>
#ifdef _OPENMP
#include <omp.h>
#endif

template <typename T>
void fRead(T *a, T *b, T *n, T *p, T input_id)
{
    std::string path_read = "/home/hexay/projects/Lab Proj/Parallel-Programing-Project/.nttdata/";
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

int main()
{
    u64 a[300000], b[300000], ab[300000];
    u64 n_, p_;

    // 使用大规模数据
    fRead<u64>(a, b, &n_, &p_, 2);

    std::ofstream outfile("performance_data.csv");
    outfile << "threads,avg_time_us,speedup,efficiency\n";

    std::cout << "开始性能测试 (n=" << n_ << ", p=" << p_ << ")" << std::endl;

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

        std::cout << "线程数: " << num_threads
                  << ", 平均时间: " << avg_time << " μs"
                  << ", 加速比: " << speedup
                  << ", 效率: " << efficiency << "%" << std::endl;
    }

    outfile.close();
    std::cout << "测试完成，数据已保存到 performance_data.csv" << std::endl;

    return 0;
}