#include "ntt/src/include/ntt.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std::chrono;

// 先定义必要的枚举
enum MulMode { MUL_NAIVE = 0, MUL_MONT = 1, MUL_BARRETT = 2 };

// 包含基础GPU函数定义（避免重复定义）
template <typename T>
__device__ T gpu_mont_reduce_correct(__uint128_t t, T mod, T neg_r_inv) {
  constexpr int word_bits = sizeof(T) * 8;
  T m = (T)t * neg_r_inv;
  __uint128_t tmp = t + (__uint128_t)m * mod;
  T res = (T)(tmp >> word_bits);
  res = res - (mod & -(res >= mod));
  return res;
}

template <typename T>
__device__ T gpu_mont_add_correct(T a_mont, T b_mont, T mod) {
  return (a_mont + b_mont >= mod) ? (a_mont + b_mont - mod) : (a_mont + b_mont);
}

template <typename T>
__device__ T gpu_mont_sub_correct(T a_mont, T b_mont, T mod) {
  return (a_mont >= b_mont) ? (a_mont - b_mont) : (a_mont + mod - b_mont);
}

template <typename T>
__device__ T gpu_mont_mul_correct(T a_mont, T b_mont, T mod, T neg_r_inv) {
  return gpu_mont_reduce_correct((__uint128_t)a_mont * b_mont, mod, neg_r_inv);
}

template <int MODE, typename T>
__device__ inline T mod_mul_unified(T a, T b, T p, T neg_r_inv, u64 mu) {
  if constexpr (MODE == MUL_MONT) {
    return gpu_mont_mul_correct(a, b, p, neg_r_inv);
  } else {
    return gpu_mont_mul_correct(a, b, p, neg_r_inv);
  }
}

template <typename T>
__global__ void pointwise_multiply_kernel_correct(T *a_mont, T *b_mont,
                                                  T *ab_mont, T n, T mod,
                                                  T neg_r_inv) {
  T idx = blockDim.x * blockIdx.x + threadIdx.x;
  if (idx < n) {
    ab_mont[idx] =
        gpu_mont_mul_correct(a_mont[idx], b_mont[idx], mod, neg_r_inv);
  }
}

template <typename T>
void generate_twiddle_table(std::vector<T> &table, T n, T p, T omega,
                            bool inverse) {
  table.resize(n);
  MontMod<T> mont(p);
  T omega_mont = mont.from_T(omega);
  if (inverse)
    omega_mont = mont.inv(omega_mont);
  T one_mont = mont.from_T(1);

  for (int level = 0; (T(1) << level) < n; ++level) {
    T mid = T(1) << level;
    T exp = (p - 1) / (mid << 1);
    T Wn_mont = mont.pow(omega_mont, exp);
    T w = one_mont;
    T offset = mid - 1;
    for (T k = 0; k < mid; ++k) {
      table[offset + k] = w;
      w = mont.mul(w, Wn_mont);
    }
  }
}

// 简化的GPU Montgomery NTT实现
template <typename T>
__global__ void ntt_stage_kernel_mont(T *a, const T *twiddles, int mid, int n,
                                      T p, T neg_r_inv) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int total = n >> 1;
  if (tid >= total)
    return;

  int step = mid << 1;
  int j = (tid / mid) * step;
  int k = tid % mid;
  int tw_idx = (mid - 1) + k;
  T w_mont = twiddles[tw_idx];

  T x = a[j + k];
  T y = gpu_mont_mul_correct(w_mont, a[j + k + mid], p, neg_r_inv);
  a[j + k] = gpu_mont_add_correct(x, y, p);
  a[j + k + mid] = gpu_mont_sub_correct(x, y, p);
}

template <typename T>
__global__ void ntt_stage_kernel_inv_mont(T *a, const T *twiddles, int mid,
                                          int n, T p, T neg_r_inv) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int total = n >> 1;
  if (tid >= total)
    return;

  int step = mid << 1;
  int j = (tid / mid) * step;
  int k = tid % mid;
  int tw_idx = (mid - 1) + k;
  T w_mont = twiddles[tw_idx];

  T x = a[j + k];
  T y = a[j + k + mid];
  a[j + k] = gpu_mont_add_correct(x, y, p);
  T temp = gpu_mont_sub_correct(x, y, p);
  a[j + k + mid] = gpu_mont_mul_correct(w_mont, temp, p, neg_r_inv);
}

// 基础GPU实现
template <typename T>
void poly_multiply_ntt_gpu_mont_basic(T *a, T *b, T *ab, T n, T p, T omega) {
  MontMod<T> montMod(p);
  T n_expanded = expand_n(2 * n - 1);

  // 计算Montgomery参数
  T neg_r_inv;
  {
    constexpr int word_bits = sizeof(T) * 8;
    T inv = 1;
    int inv_iterations = (word_bits <= 64) ? 6 : 7;
    for (int i = 0; i < inv_iterations; ++i) {
      inv = (T)((__uint128_t)inv * (2 - (__uint128_t)p * inv));
    }
    neg_r_inv = -inv;
  }

  // GPU内存分配
  T *a_mont_gpu, *b_mont_gpu, *ab_mont_gpu, *tw_gpu, *tw_inv_gpu;
  cudaMalloc(&a_mont_gpu, n_expanded * sizeof(T));
  cudaMalloc(&b_mont_gpu, n_expanded * sizeof(T));
  cudaMalloc(&ab_mont_gpu, n_expanded * sizeof(T));
  cudaMalloc(&tw_gpu, n_expanded * sizeof(T));
  cudaMalloc(&tw_inv_gpu, n_expanded * sizeof(T));

  // 数据预处理
  T *a_expanded = expand_a(a, n, n_expanded);
  T *b_expanded = expand_a(b, n, n_expanded);
  bit_reverse_permute(a_expanded, n_expanded);
  bit_reverse_permute(b_expanded, n_expanded);

  // Montgomery转换
  T *a_mont_host = new T[n_expanded];
  T *b_mont_host = new T[n_expanded];
  for (T i = 0; i < n_expanded; i++) {
    a_mont_host[i] = montMod.from_T(a_expanded[i]);
    b_mont_host[i] = montMod.from_T(b_expanded[i]);
  }

  // 数据传输
  cudaMemcpy(a_mont_gpu, a_mont_host, n_expanded * sizeof(T),
             cudaMemcpyHostToDevice);
  cudaMemcpy(b_mont_gpu, b_mont_host, n_expanded * sizeof(T),
             cudaMemcpyHostToDevice);

  // 旋转因子
  std::vector<T> tw_host(n_expanded), tw_inv_host(n_expanded);
  generate_twiddle_table(tw_host, n_expanded, p, omega, false);
  generate_twiddle_table(tw_inv_host, n_expanded, p, omega, true);
  cudaMemcpy(tw_gpu, tw_host.data(), n_expanded * sizeof(T),
             cudaMemcpyHostToDevice);
  cudaMemcpy(tw_inv_gpu, tw_inv_host.data(), n_expanded * sizeof(T),
             cudaMemcpyHostToDevice);

  // NTT前向变换
  const int BLOCK = 256;
  for (int level = 0; (T(1) << level) < n_expanded; ++level) {
    T mid = T(1) << level;
    int total_butterfly = n_expanded >> 1;
    int GRID = (total_butterfly + BLOCK - 1) / BLOCK;

    ntt_stage_kernel_mont<<<GRID, BLOCK>>>(a_mont_gpu, tw_gpu, mid, n_expanded,
                                           p, neg_r_inv);
    ntt_stage_kernel_mont<<<GRID, BLOCK>>>(b_mont_gpu, tw_gpu, mid, n_expanded,
                                           p, neg_r_inv);
  }
  cudaDeviceSynchronize();

  // 点乘
  dim3 blockSize(256);
  dim3 gridSize((n_expanded + blockSize.x - 1) / blockSize.x);
  pointwise_multiply_kernel_correct<<<gridSize, blockSize>>>(
      a_mont_gpu, b_mont_gpu, ab_mont_gpu, n_expanded, p, neg_r_inv);
  cudaDeviceSynchronize();

  // NTT逆变换
  for (int level = (int)std::log2(n_expanded) - 1; level >= 0; --level) {
    T mid = T(1) << level;
    int total_butterfly = n_expanded >> 1;
    int GRID = (total_butterfly + BLOCK - 1) / BLOCK;
    ntt_stage_kernel_inv_mont<<<GRID, BLOCK>>>(ab_mont_gpu, tw_inv_gpu, mid,
                                               n_expanded, p, neg_r_inv);
  }
  cudaDeviceSynchronize();

  // 结果处理
  T *temp_result = new T[n_expanded];
  cudaMemcpy(temp_result, ab_mont_gpu, n_expanded * sizeof(T),
             cudaMemcpyDeviceToHost);

  T inv_n_mont = montMod.inv(montMod.from_T(n_expanded));
  for (T i = 0; i < n_expanded; i++) {
    temp_result[i] = montMod.mul(temp_result[i], inv_n_mont);
    temp_result[i] = montMod.to_T(temp_result[i]);
  }

  bit_reverse_permute(temp_result, n_expanded);
  for (T i = 0; i < 2 * n - 1; i++) {
    ab[i] = temp_result[i];
  }

  // 清理
  delete[] a_expanded;
  delete[] b_expanded;
  delete[] a_mont_host;
  delete[] b_mont_host;
  delete[] temp_result;
  cudaFree(a_mont_gpu);
  cudaFree(b_mont_gpu);
  cudaFree(ab_mont_gpu);
  cudaFree(tw_gpu);
  cudaFree(tw_inv_gpu);
}

// 现在包含优化版本（已修复重复定义）
#include "ntt/src/include/CUDA/ntt_gpu_optimized.cu"

template <typename T> void generate_test_data(T *a, T *b, T n) {
  for (T i = 0; i < n; i++) {
    a[i] = (i * 17 + 7) % 1000 + 1;
    b[i] = (i * 23 + 11) % 1000 + 1;
  }
}

void run_optimization_test() {
  const uint32_t p = 998244353;
  const uint32_t omega = 3;

  std::vector<uint32_t> test_sizes = {4096, 16384, 65536, 131072};

  std::cout << "=== GPU优化效果测试 ===" << std::endl;
  std::cout << std::setw(10) << "Size" << std::setw(15) << "CPU时间(μs)"
            << std::setw(15) << "GPU基础(μs)" << std::setw(15) << "GPU优化(μs)"
            << std::setw(12) << "基础加速比" << std::setw(12) << "优化加速比"
            << std::setw(12) << "改进倍数" << std::endl;
  std::cout << std::string(90, '-') << std::endl;

  for (uint32_t n : test_sizes) {
    uint32_t *a = new uint32_t[n];
    uint32_t *b = new uint32_t[n];
    uint32_t *cpu_result = new uint32_t[2 * n - 1]();
    uint32_t *gpu_basic_result = new uint32_t[2 * n - 1]();
    uint32_t *gpu_opt_result = new uint32_t[2 * n - 1]();

    generate_test_data(a, b, n);

    // CPU基线测试
    uint32_t *a_copy = new uint32_t[n];
    uint32_t *b_copy = new uint32_t[n];
    std::copy(a, a + n, a_copy);
    std::copy(b, b + n, b_copy);

    auto cpu_start = high_resolution_clock::now();
    poly_multiply_ntt(a_copy, b_copy, cpu_result, n, p, omega);
    auto cpu_end = high_resolution_clock::now();
    double cpu_time = duration_cast<microseconds>(cpu_end - cpu_start).count();

    delete[] a_copy;
    delete[] b_copy;

    // GPU基础版本测试
    a_copy = new uint32_t[n];
    b_copy = new uint32_t[n];
    std::copy(a, a + n, a_copy);
    std::copy(b, b + n, b_copy);

    auto gpu_basic_start = high_resolution_clock::now();
    poly_multiply_ntt_gpu_mont_basic(a_copy, b_copy, gpu_basic_result, n, p,
                                     omega);
    auto gpu_basic_end = high_resolution_clock::now();
    double gpu_basic_time =
        duration_cast<microseconds>(gpu_basic_end - gpu_basic_start).count();

    delete[] a_copy;
    delete[] b_copy;

    // GPU优化版本测试
    a_copy = new uint32_t[n];
    b_copy = new uint32_t[n];
    std::copy(a, a + n, a_copy);
    std::copy(b, b + n, b_copy);

    auto gpu_opt_start = high_resolution_clock::now();
    poly_multiply_ntt_gpu_mont_optimized(a_copy, b_copy, gpu_opt_result, n, p,
                                         omega);
    auto gpu_opt_end = high_resolution_clock::now();
    double gpu_opt_time =
        duration_cast<microseconds>(gpu_opt_end - gpu_opt_start).count();

    delete[] a_copy;
    delete[] b_copy;

    // 验证正确性
    bool basic_correct = true, opt_correct = true;
    for (uint32_t i = 0; i < 2 * n - 1; i++) {
      if (cpu_result[i] != gpu_basic_result[i])
        basic_correct = false;
      if (cpu_result[i] != gpu_opt_result[i])
        opt_correct = false;
    }

    // 计算性能指标
    double basic_speedup = cpu_time / gpu_basic_time;
    double opt_speedup = cpu_time / gpu_opt_time;
    double improvement = gpu_basic_time / gpu_opt_time;

    std::cout << std::setw(10) << n << std::setw(15) << std::fixed
              << std::setprecision(1) << cpu_time << std::setw(15) << std::fixed
              << std::setprecision(1) << gpu_basic_time << std::setw(15)
              << std::fixed << std::setprecision(1) << gpu_opt_time
              << std::setw(12) << std::fixed << std::setprecision(2)
              << basic_speedup << std::setw(12) << std::fixed
              << std::setprecision(2) << opt_speedup << std::setw(12)
              << std::fixed << std::setprecision(2) << improvement << " ["
              << (basic_correct ? "✓" : "✗") << "," << (opt_correct ? "✓" : "✗")
              << "]" << std::endl;

    delete[] a;
    delete[] b;
    delete[] cpu_result;
    delete[] gpu_basic_result;
    delete[] gpu_opt_result;
  }
}

int main() {
  // 检查CUDA设备
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  if (deviceCount == 0) {
    std::cerr << "没有找到CUDA设备!" << std::endl;
    return 1;
  }

  // 显示GPU信息
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, 0);
  std::cout << "GPU设备: " << prop.name << std::endl;
  std::cout << "SM数量: " << prop.multiProcessorCount << std::endl;
  std::cout << std::endl;

  try {
    run_optimization_test();
  } catch (const std::exception &e) {
    std::cerr << "测试错误: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}