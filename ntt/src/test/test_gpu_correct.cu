#include "src/include/CUDA/ntt.h"
#include <iostream>
#include <chrono>

void test_small_case() {
    u64 a[] = {1, 2, 3, 4};
    u64 b[] = {5, 6, 7, 8};
    u64 ab[7] = {0};
    u64 n = 4;
    u64 p = 7340033;
    u64 omega = 3;
    
    std::cout << "Testing small case: a=[1,2,3,4], b=[5,6,7,8], p=" << p << std::endl;
    
    auto start = std::chrono::high_resolution_clock::now();
    poly_multiply_ntt_gpu_correct(a, b, ab, n, p, omega);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "GPU Correct version time: " << duration.count() << " μs" << std::endl;
    
    std::cout << "Result: [";
    for (int i = 0; i < 7; i++) {
        std::cout << ab[i];
        if (i < 6) std::cout << ", ";
    }
    std::cout << "]" << std::endl;
    
    // 期望结果：[5, 16, 34, 60, 61, 52, 32]
    u64 expected[] = {5, 16, 34, 60, 61, 52, 32};
    bool correct = true;
    for (int i = 0; i < 7; i++) {
        if (ab[i] != expected[i]) {
            correct = false;
            break;
        }
    }
    
    if (correct) {
        std::cout << "✅ CORRECT! Matches expected result." << std::endl;
    } else {
        std::cout << "❌ INCORRECT! Expected: [5, 16, 34, 60, 61, 52, 32]" << std::endl;
    }
}

int main() {
    std::cout << "=== Testing GPU Correct NTT Implementation ===" << std::endl;
    
    // 检查GPU设备
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    std::cout << "Found " << deviceCount << " CUDA device(s)" << std::endl;
    
    if (deviceCount == 0) {
        std::cout << "No CUDA devices found!" << std::endl;
        return 1;
    }
    
    // 设备信息
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    std::cout << "Using device: " << prop.name << std::endl;
    std::cout << "Memory: " << prop.totalGlobalMem / (1024*1024) << " MB" << std::endl;
    
    // 运行测试
    test_small_case();
    
    return 0;
} 