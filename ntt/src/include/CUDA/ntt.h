#pragma once

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
template <typename T>
__device__ T gpu_mont_add(T a_mont, T b_mont, T mod);

template <typename T>
__device__ T gpu_mont_sub(T a_mont, T b_mont, T mod);

template <typename T>
__device__ T gpu_mont_mul(T a_mont, T b_mont, T mod, T neg_r_inv);

template <typename T>
__device__ T gpu_mont_reduce(uint64_t t, T mod, T neg_r_inv);

// NTT核函数声明
template <typename T>
__global__ void ntt_forward_kernel_basic(T *a_mont, T n, T p, T *twiddle_factors, int mid_level);

template <typename T>
__global__ void ntt_forward_kernel_optimized(T *a_mont, T n, T p, T *twiddle_factors, int mid_level);

template <typename T>
__global__ void ntt_inverse_kernel_basic(T *a_mont, T n, T p, T *twiddle_factors, int mid_level);

template <typename T>
__global__ void ntt_inverse_kernel_optimized(T *a_mont, T n, T p, T *twiddle_factors, int mid_level);

template <typename T>
__global__ void pointwise_multiply_kernel(T *a_mont, T *b_mont, T *ab_mont, T n, T mod, T neg_r_inv);

// GPU版本的多项式乘法函数
template <typename T>
void poly_multiply_ntt_gpu_basic(T *a, T *b, T *ab, T n, T p, T omega = 3);

template <typename T>
void poly_multiply_ntt_gpu_optimized(T *a, T *b, T *ab, T n, T p, T omega = 3);

template <typename T>
void poly_multiply_ntt_gpu(T *a, T *b, T *ab, T n, T p, T omega = 3);

// GPU 多模乘版本接口
template <typename T>
void poly_multiply_ntt_gpu_naive(T *a, T *b, T *ab, T n, T p, T omega = 3);

template <typename T>
void poly_multiply_ntt_gpu_mont(T *a, T *b, T *ab, T n, T p, T omega = 3);

template <typename T>
void poly_multiply_ntt_gpu_barrett(T *a, T *b, T *ab, T n, T p, T omega = 3);

template <typename T>
void copy_to_gpu_and_expand(T *host_data, T **gpu_data, T n, T n_expanded);

template <typename T>
void copy_from_gpu_and_shrink(T *gpu_data, T *host_data, T n_expanded, T n);

// 辅助函数
template <typename T>
void precompute_twiddle_factors_gpu(T **twiddle_factors_gpu, T n, T p, T omega_mont, bool is_inverse = false);

// 模板函数的实现将在编译时包含