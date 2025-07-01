#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cstdlib>
#include <cassert>
#include <fstream>

#include "src/include/ntt.h"

// 确保GPU版本的声明可见
template <typename T>
void poly_multiply_ntt_gpu_naive(T *a, T *b, T *ab, T n, T p, T omega);

template <typename T>
void poly_multiply_ntt_gpu_mont(T *a, T *b, T *ab, T n, T p, T omega);

template <typename T>
void poly_multiply_ntt_gpu_barrett(T *a, T *b, T *ab, T n, T p, T omega);

template <typename T>
void poly_multiply_ntt_gpu(T *a, T *b, T *ab, T n, T p, T omega);

// CUDA工具
#define CHECK_CUDA_CORRECT(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA error: " << cudaGetErrorString(err) << std::endl; \
            exit(1); \
        } \
    } while(0)

void cuda_check_and_reset() {
    CHECK_CUDA_CORRECT(cudaDeviceSynchronize());
    CHECK_CUDA_CORRECT(cudaDeviceReset());
}

// 测试用例生成
template <typename T>
void generate_random_poly(T *a, T n, T max_val) {
    for (T i = 0; i < n; i++) {
        a[i] = rand() % max_val;
    }
}

// 结果验证
template <typename T>
bool verify_results(const T *expected, const T *actual, T len, const std::string &test_name) {
    for (T i = 0; i < len; i++) {
        if (expected[i] != actual[i]) {
            std::cout << "❌ " << test_name << " FAILED at index " << i 
                      << ": expected=" << expected[i] << ", actual=" << actual[i] << std::endl;
            return false;
        }
    }
    std::cout << "✅ " << test_name << " PASSED" << std::endl;
    return true;
}

// 高精度计时
template <typename Func>
double measure_time_us(Func func, int repeat = 5) {
    double total_time = 0.0;
    for (int i = 0; i < repeat; i++) {
        auto start = std::chrono::high_resolution_clock::now();
        func();
        CHECK_CUDA_CORRECT(cudaDeviceSynchronize()); // 确保GPU完成
        auto end = std::chrono::high_resolution_clock::now();
        
        double elapsed_us = std::chrono::duration<double, std::micro>(end - start).count();
        total_time += elapsed_us;
    }
    return total_time / repeat;
}

// 性能对比结构
struct PerfResult {
    std::string name;
    double time_us;
    bool correct;
    double speedup;
};

// 主测试函数
template <typename T>
void run_performance_comparison(T n, T p, T omega, const std::string &test_name) {
    std::cout << "\n=== " << test_name << " (n=" << n << ") ===" << std::endl;
    
    // 生成随机测试数据
    T *a = new T[n];
    T *b = new T[n];
    generate_random_poly(a, n, p / 100);
    generate_random_poly(b, n, p / 100);
    
    // 结果数组
    T result_len = 2 * n - 1;
    T *cpu_result = new T[result_len]();
    T *gpu_naive_result = new T[result_len]();
    T *gpu_mont_result = new T[result_len]();
    T *gpu_barrett_result = new T[result_len]();
    
    std::vector<PerfResult> results;
    
    // === CPU基准测试 ===
    auto cpu_func = [&]() {
        std::fill(cpu_result, cpu_result + result_len, 0);
        poly_multiply_ntt(a, b, cpu_result, n, p, omega);
    };
    double cpu_time = measure_time_us(cpu_func);
    results.push_back({"CPU", cpu_time, true, 1.0});
    
    // === GPU朴素模乘测试 ===
    bool naive_correct = false;
    double naive_time = 0.0;
    try {
        auto naive_func = [&]() {
            std::fill(gpu_naive_result, gpu_naive_result + result_len, 0);
            poly_multiply_ntt_gpu_naive(a, b, gpu_naive_result, n, p, omega);
        };
        naive_time = measure_time_us(naive_func);
        naive_correct = verify_results(cpu_result, gpu_naive_result, result_len, "GPU Naive");
    } catch (...) {
        std::cout << "❌ GPU Naive 运行失败" << std::endl;
    }
    results.push_back({"GPU_Naive", naive_time, naive_correct, cpu_time / naive_time});
    
    // === GPU Montgomery模乘测试 ===
    bool mont_correct = false;
    double mont_time = 0.0;
    try {
        auto mont_func = [&]() {
            std::fill(gpu_mont_result, gpu_mont_result + result_len, 0);
            poly_multiply_ntt_gpu_mont(a, b, gpu_mont_result, n, p, omega);
        };
        mont_time = measure_time_us(mont_func);
        mont_correct = verify_results(cpu_result, gpu_mont_result, result_len, "GPU Montgomery");
    } catch (...) {
        std::cout << "❌ GPU Montgomery 运行失败" << std::endl;
    }
    results.push_back({"GPU_Montgomery", mont_time, mont_correct, cpu_time / mont_time});
    
    // === GPU Barrett模乘测试 ===
    bool barrett_correct = false;
    double barrett_time = 0.0;
    try {
        auto barrett_func = [&]() {
            std::fill(gpu_barrett_result, gpu_barrett_result + result_len, 0);
            poly_multiply_ntt_gpu_barrett(a, b, gpu_barrett_result, n, p, omega);
        };
        barrett_time = measure_time_us(barrett_func);
        barrett_correct = verify_results(cpu_result, gpu_barrett_result, result_len, "GPU Barrett");
    } catch (...) {
        std::cout << "❌ GPU Barrett 运行失败" << std::endl;
    }
    results.push_back({"GPU_Barrett", barrett_time, barrett_correct, cpu_time / barrett_time});
    
    // === 结果统计 ===
    std::cout << "\n性能对比结果:" << std::endl;
    std::cout << std::left << std::setw(15) << "版本" 
              << std::setw(12) << "时间(μs)" 
              << std::setw(10) << "正确性"
              << std::setw(10) << "相对速度" << std::endl;
    std::cout << std::string(50, '-') << std::endl;
    
    for (const auto &result : results) {
        std::cout << std::left << std::setw(15) << result.name
                  << std::setw(12) << std::fixed << std::setprecision(1) << result.time_us
                  << std::setw(10) << (result.correct ? "✅" : "❌")
                  << std::setw(10) << std::fixed << std::setprecision(2) << result.speedup << "x"
                  << std::endl;
    }
    
    // 保存CSV数据
    static bool csv_header_written = false;
    std::ofstream csv_file("perf_results.csv", csv_header_written ? std::ios::app : std::ios::out);
    if (!csv_header_written) {
        csv_file << "TestName,n,p,CPU_time,GPU_Naive_time,GPU_Mont_time,GPU_Barrett_time,";
        csv_file << "GPU_Naive_speedup,GPU_Mont_speedup,GPU_Barrett_speedup,";
        csv_file << "GPU_Naive_correct,GPU_Mont_correct,GPU_Barrett_correct" << std::endl;
        csv_header_written = true;
    }
    
    csv_file << test_name << "," << n << "," << p << ",";
    csv_file << cpu_time << "," << naive_time << "," << mont_time << "," << barrett_time << ",";
    csv_file << (cpu_time / naive_time) << "," << (cpu_time / mont_time) << "," << (cpu_time / barrett_time) << ",";
    csv_file << naive_correct << "," << mont_correct << "," << barrett_correct << std::endl;
    
    // 清理内存
    delete[] a;
    delete[] b;
    delete[] cpu_result;
    delete[] gpu_naive_result;
    delete[] gpu_mont_result;
    delete[] gpu_barrett_result;
}

int main() {
    srand(42); // 确保可重现性
    
    std::cout << "GPU NTT 模乘算法性能对比实验" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // 检查CUDA环境
    int deviceCount;
    CHECK_CUDA_CORRECT(cudaGetDeviceCount(&deviceCount));
    std::cout << "发现 " << deviceCount << " 个CUDA设备" << std::endl;
    
    if (deviceCount == 0) {
        std::cout << "未找到CUDA设备，退出。" << std::endl;
        return 1;
    }
    
    // 删除旧的CSV文件
    std::remove("perf_results.csv");
    
    // === 32位测试用例 ===
    std::cout << "\n=== 32位整数测试 ===" << std::endl;
    run_performance_comparison<u32>(4, 7340033, 5, "Small_u32");
    run_performance_comparison<u32>(256, 7340033, 5, "Medium_u32");
    run_performance_comparison<u32>(4096, 7340033, 5, "Large_u32");
    run_performance_comparison<u32>(65536, 7340033, 5, "XLarge_u32");
    
    // === 64位测试用例 ===
    std::cout << "\n=== 64位整数测试 ===" << std::endl;
    u64 p64 = 2305843009213693951ULL; // 2^61 - 1
    u64 omega64 = 1753635133440165772ULL;
    
    run_performance_comparison<u64>(4, p64, omega64, "Small_u64");
    run_performance_comparison<u64>(256, p64, omega64, "Medium_u64");
    run_performance_comparison<u64>(4096, p64, omega64, "Large_u64");
    run_performance_comparison<u64>(65536, p64, omega64, "XLarge_u64");
    
    std::cout << "\n实验完成！结果已保存到 perf_results.csv" << std::endl;
    
    cuda_check_and_reset();
    return 0;
} 