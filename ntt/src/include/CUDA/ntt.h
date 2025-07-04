#pragma once

#include <vector>
#include "../general/utils.h"
#include "../general/op.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// CUDA错误检查宏
#define CHECK_CUDA(call)                                                                               \
    do                                                                                                 \
    {                                                                                                  \
        cudaError_t err = call;                                                                        \
        if (err != cudaSuccess)                                                                        \
        {                                                                                              \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
            exit(1);                                                                                   \
        }                                                                                              \
    } while (0)

// GPU设备函数声明
__device__ u32 gpu_mont_add(u32 a_mont, u32 b_mont, u32 mod);

__device__ u32 gpu_mont_sub(u32 a_mont, u32 b_mont, u32 mod);

__device__ u32 gpu_mont_mul(u32 a_mont, u32 b_mont, u32 mod, u32 neg_r_inv);

__device__ u32 gpu_mont_reduce(uint64_t t, u32 mod, u32 neg_r_inv);

// NTT核函数声明
__global__ void ntt_forward_kernel_basic(u32 *a_mont, u32 n, u32 p, u32 *twiddle_factors, int mid_level);

__global__ void ntt_forward_kernel_optimized(u32 *a_mont, u32 n, u32 p, u32 *twiddle_factors, int mid_level);

__global__ void ntt_inverse_kernel_basic(u32 *a_mont, u32 n, u32 p, u32 *twiddle_factors, int mid_level);

__global__ void ntt_inverse_kernel_optimized(u32 *a_mont, u32 n, u32 p, u32 *twiddle_factors, int mid_level);

__global__ void pointwise_multiply_kernel(u32 *a_mont, u32 *b_mont, u32 *ab_mont, u32 n, u32 mod, u32 neg_r_inv);

// GPU版本的多项式乘法函数
void poly_multiply_ntt_gpu_basic(u32 *a, u32 *b, u32 *ab, u32 n, u32 p, u32 omega = 3);

void poly_multiply_ntt_gpu(u32 *a, u32 *b, u32 *ab, u32 n, u32 p, u32 omega = 3);

// GPU 多模乘版本接口
void poly_multiply_ntt_gpu_naive(u32 *a, u32 *b, u32 *ab, u32 n, u32 p, u32 omega = 3);

void poly_multiply_ntt_gpu_mont(u32 *a, u32 *b, u32 *ab, u32 n, u32 p, u32 omega = 3);

void poly_multiply_ntt_gpu_barrett(u32 *a, u32 *b, u32 *ab, u32 n, u32 p, u32 omega = 3);

// 优化版本接口 - Lab5.tex优化策略实现
void poly_multiply_ntt_gpu_mont_optimized(u32 *a, u32 *b, u32 *ab, u32 n, u32 p, u32 omega = 3);

void poly_multiply_ntt_gpu_naive_optimized(u32 *a, u32 *b, u32 *ab, u32 n, u32 p, u32 omega = 3);

void poly_multiply_ntt_gpu_barrett_optimized(u32 *a, u32 *b, u32 *ab, u32 n, u32 p, u32 omega = 3);

// 旋转因子预计算函数声明 (u32版本)
void generate_twiddle_table_u32(std::vector<u32> &table, u32 n, u32 p, u32 omega, bool inverse);

void copy_to_gpu_and_expand(u32 *host_data, u32 **gpu_data, u32 n, u32 n_expanded);

void copy_from_gpu_and_shrink(u32 *gpu_data, u32 *host_data, u32 n_expanded, u32 n);

// 辅助函数
void precompute_twiddle_factors_gpu(u32 **twiddle_factors_gpu, u32 n, u32 p, u32 omega_mont, bool is_inverse = false);

// 模板函数的实现将在编译时包含