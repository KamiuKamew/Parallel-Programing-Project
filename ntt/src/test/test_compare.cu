#include "src/include/ntt.h"
#include "src/include/CUDA/ntt_gpu.cu"
#include <iostream>
#include <vector>
#include <cstring>

using namespace std;

int main() {
    cout << "Comparing CPU vs GPU NTT implementation..." << endl;
    
    // 简单的测试数据
    const u64 n = 4;
    const u64 p = 7340033;  // 第一个测试模数
    const u64 omega = 3;    // 原根
    
    u64 a[n] = {1, 2, 3, 4};
    u64 b[n] = {5, 6, 7, 8};
    u64 ab_cpu[2*n-1] = {0};
    u64 ab_gpu[2*n-1] = {0};
    
    cout << "Input a: ";
    for (int i = 0; i < n; i++) {
        cout << a[i] << " ";
    }
    cout << endl;
    
    cout << "Input b: ";
    for (int i = 0; i < n; i++) {
        cout << b[i] << " ";
    }
    cout << endl;
    
    // CPU版本
    try {
        poly_multiply_ntt(a, b, ab_cpu, n, p, omega);
        cout << "CPU result: ";
        for (int i = 0; i < 2*n-1; i++) {
            cout << ab_cpu[i] << " ";
        }
        cout << endl;
    } catch (const exception& e) {
        cout << "CPU Error: " << e.what() << endl;
    }
    
    // GPU版本
    try {
        poly_multiply_ntt_gpu_basic(a, b, ab_gpu, n, p, omega);
        cout << "GPU result: ";
        for (int i = 0; i < 2*n-1; i++) {
            cout << ab_gpu[i] << " ";
        }
        cout << endl;
    } catch (const exception& e) {
        cout << "GPU Error: " << e.what() << endl;
    }
    
    // 朴素算法作为参考
    u64 ab_naive[2*n-1] = {0};
    for (u64 i = 0; i < n; ++i) {
        for (u64 j = 0; j < n; ++j) {
            ab_naive[i + j] = (1LL * a[i] * b[j] % p + ab_naive[i + j]) % p;
        }
    }
    
    cout << "Naive result: ";
    for (int i = 0; i < 2*n-1; i++) {
        cout << ab_naive[i] << " ";
    }
    cout << endl;
    
    // 对比结果
    bool cpu_correct = true, gpu_correct = true;
    for (int i = 0; i < 2*n-1; i++) {
        if (ab_cpu[i] != ab_naive[i]) cpu_correct = false;
        if (ab_gpu[i] != ab_naive[i]) gpu_correct = false;
    }
    
    cout << "CPU correctness: " << (cpu_correct ? "PASS" : "FAIL") << endl;
    cout << "GPU correctness: " << (gpu_correct ? "PASS" : "FAIL") << endl;
    
    return 0;
} 