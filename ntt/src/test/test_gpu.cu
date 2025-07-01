#include "src/include/CUDA/ntt_gpu.cu"
#include <iostream>
#include <vector>

using namespace std;

int main() {
    cout << "Testing GPU NTT implementation..." << endl;
    
    // 简单的测试数据
    const u64 n = 4;
    const u64 p = 7340033;  // 第一个测试模数
    const u64 omega = 3;    // 原根
    
    u64 a[n] = {1, 2, 3, 4};
    u64 b[n] = {5, 6, 7, 8};
    u64 ab[2*n-1] = {0};
    
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
    
    try {
        // 调用GPU版本的多项式乘法
        poly_multiply_ntt_gpu_basic(a, b, ab, n, p, omega);
        
        cout << "GPU result: ";
        for (int i = 0; i < 2*n-1; i++) {
            cout << ab[i] << " ";
        }
        cout << endl;
        
        cout << "GPU NTT test completed successfully!" << endl;
        
    } catch (const exception& e) {
        cout << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
} 