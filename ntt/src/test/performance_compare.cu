#include "src/include/CUDA/ntt.h" 
#include "src/include/ntt.h"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <cstring>

void compare_performance() {
    std::cout << "=== CPU vs GPU Performance Comparison ===" << std::endl;
    std::cout << std::setw(8) << "n" << std::setw(15) << "p" 
              << std::setw(12) << "CPU(μs)" << std::setw(12) << "GPU(μs)" 
              << std::setw(10) << "Speedup" << std::setw(10) << "Status" << std::endl;
    std::cout << std::string(70, '-') << std::endl;
    
    // 小规模测试
    {
        u64 a[] = {1, 2, 3, 4};
        u64 b[] = {5, 6, 7, 8};
        u64 ab_cpu[7] = {0}, ab_gpu[7] = {0};
        u64 n = 4, p = 7340033, omega = 3;
        
        // CPU测试
        auto start = std::chrono::high_resolution_clock::now();
        poly_multiply_ntt<u64>(a, b, ab_cpu, n, p, omega);
        auto end = std::chrono::high_resolution_clock::now();
        auto cpu_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        
        // GPU测试
        start = std::chrono::high_resolution_clock::now();
        poly_multiply_ntt_gpu_correct<u64>(a, b, ab_gpu, n, p, omega);
        end = std::chrono::high_resolution_clock::now();
        auto gpu_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        
        // 验证正确性
        bool correct = true;
        for (int i = 0; i < 7; i++) {
            if (ab_cpu[i] != ab_gpu[i]) {
                correct = false;
                break;
            }
        }
        
        double speedup = (double)cpu_time / gpu_time;
        std::cout << std::setw(8) << n << std::setw(15) << p 
                  << std::setw(12) << cpu_time << std::setw(12) << gpu_time 
                  << std::setw(10) << std::fixed << std::setprecision(2) << speedup 
                  << std::setw(10) << (correct ? "✅" : "❌") << std::endl;
    }
    
    // 从实际测试文件读取大规模数据
    std::cout << "\n--- Reading test data from files ---" << std::endl;
    
    for (int test_id = 1; test_id <= 4; test_id++) {
        std::string data_path = "/home/hexay/projects/Lab Proj/Parallel-Programing-Project/.nttdata/" 
                               + std::to_string(test_id) + ".in";
        
        std::ifstream fin(data_path);
        if (!fin.is_open()) {
            std::cout << "Cannot open " << data_path << std::endl;
            continue;
        }
        
        u64 n, p;
        fin >> n >> p;
        
        u64 *a = new u64[n];
        u64 *b = new u64[n];
        u64 *ab_cpu = new u64[2*n-1];
        u64 *ab_gpu = new u64[2*n-1];
        
        for (u64 i = 0; i < n; i++) fin >> a[i];
        for (u64 i = 0; i < n; i++) fin >> b[i];
        fin.close();
        
        memset(ab_cpu, 0, (2*n-1) * sizeof(u64));
        memset(ab_gpu, 0, (2*n-1) * sizeof(u64));
        
        // CPU测试
        auto start = std::chrono::high_resolution_clock::now();
        poly_multiply_ntt<u64>(a, b, ab_cpu, n, p, 3);
        auto end = std::chrono::high_resolution_clock::now();
        auto cpu_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        
        // GPU测试  
        start = std::chrono::high_resolution_clock::now();
        poly_multiply_ntt_gpu_correct<u64>(a, b, ab_gpu, n, p, 3);
        end = std::chrono::high_resolution_clock::now();
        auto gpu_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        
        // 验证正确性
        bool correct = true;
        for (u64 i = 0; i < 2*n-1; i++) {
            if (ab_cpu[i] != ab_gpu[i]) {
                correct = false;
                break;
            }
        }
        
        double speedup = (double)cpu_time / gpu_time;
        std::cout << std::setw(8) << n << std::setw(15) << p 
                  << std::setw(12) << cpu_time << std::setw(12) << gpu_time 
                  << std::setw(10) << speedup 
                  << std::setw(10) << (correct ? "✅" : "❌") << std::endl;
        
        delete[] a;
        delete[] b; 
        delete[] ab_cpu;
        delete[] ab_gpu;
    }
}

int main() {
    compare_performance();
    return 0;
} 