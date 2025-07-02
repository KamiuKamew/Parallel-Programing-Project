#pragma once

#include "../general/op.h"
#include "../general/type.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// GPU常量内存用于存储旋转因子 - 调整大小避免超出限制
__constant__ u64 gpu_twiddle_factors[4096]; // 支持最大n=4096的旋转因子

// GPU错误检查宏
#define CUDA_CHECK(call)                                                                       \
    do                                                                                         \
    {                                                                                          \
        cudaError_t err = call;                                                                \
        if (err != cudaSuccess)                                                                \
        {                                                                                      \
            printf("CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
            exit(1);                                                                           \
        }                                                                                      \
    } while (0)

// GPU设备函数声明
__device__ u64 gpu_mod_add(u64 a, u64 b, u64 mod);
__device__ u64 gpu_mod_sub(u64 a, u64 b, u64 mod);
__device__ u64 gpu_mod_mul_naive(u64 a, u64 b, u64 mod);
__device__ u64 gpu_mod_pow(u64 base, u64 exp, u64 mod);

// GPU kernel函数声明
__global__ void ntt_forward_gpu_kernel(u64 *a, u64 n, u64 mid, u64 p, u64 Wn);
__global__ void ntt_inverse_gpu_kernel(u64 *a, u64 n, u64 mid, u64 p, u64 Wn);

// 主机端接口函数声明
void ntt_forward_gpu(u64 *a, u64 n, u64 p, u64 omega);
void ntt_inverse_gpu(u64 *a, u64 n, u64 p, u64 omega);
void poly_multiply_ntt_gpu(u64 *a, u64 *b, u64 *result, u64 n, u64 p, u64 omega);

// 旋转因子预处理函数
void precompute_twiddle_factors_gpu(u64 n, u64 p, u64 omega);