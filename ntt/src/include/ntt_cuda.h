#pragma once

#include "general/utils.h"
#include "transform.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

/**
 * @brief 使用CUDA的NTT优化多项式乘法 - 第一步：最简化版本
 *
 * 这个版本使用最简单的GPU朴素乘法，主要用于验证CUDA环境能正常工作
 *
 * @param a 多项式系数
 * @param b 多项式系数
 * @param ab 结果
 * @param n 多项式长度
 * @param p 模数（质数）
 */
void poly_multiply_ntt_cuda_simple(u64 *a, u64 *b, u64 *ab, u64 n, u64 p,
                                   u64 OMEGA = 3);

/**
 * @brief 使用CUDA的NTT优化多项式乘法 - 第一步：基础串行版本
 *
 * 这个版本使用CUDA内存管理但保持串行计算，主要用于验证CUDA环境和内存操作的正确性
 *
 * @param a 多项式系数
 * @param b 多项式系数
 * @param ab 结果
 * @param n 多项式长度
 * @param p 模数（质数）
 */
void poly_multiply_ntt_cuda_serial(u64 *a, u64 *b, u64 *ab, u64 n, u64 p,
                                   u64 OMEGA = 3);

/**
 * @brief 使用CUDA的NTT优化多项式乘法 - 第二步：便于并行化的串行版本
 *
 * 这个版本重构了循环结构，为GPU并行化做准备，但仍在GPU上串行执行
 * 预处理旋转因子，调整数据布局便于GPU访问
 *
 * @param a 多项式系数
 * @param b 多项式系数
 * @param ab 结果
 * @param n 多项式长度
 * @param p 模数（质数）
 */
void poly_multiply_ntt_cuda_parallel_ready(u64 *a, u64 *b, u64 *ab, u64 n,
                                           u64 p, u64 OMEGA = 3);

/**
 * @brief 使用CUDA的NTT优化多项式乘法 - 第三步：完整GPU并行版本
 *
 * 这个版本实现了完整的GPU并行化，将串行的蝶形运算并行化到GPU线程
 * 每个线程处理一个独立的蝶形运算，真正发挥GPU并行计算能力
 *
 * @param a 多项式系数
 * @param b 多项式系数
 * @param ab 结果
 * @param n 多项式长度
 * @param p 模数（质数）
 */
void poly_multiply_ntt_cuda_parallel(u64 *a, u64 *b, u64 *ab, u64 n, u64 p,
                                     u64 OMEGA = 3);

/**
 * @brief 使用CUDA的NTT优化多项式乘法 - 完整GPU并行版本（模板）
 *
 * 这个版本实现了完整的GPU并行化，对第二层和第三层循环进行并行化
 *
 * @param a 多项式系数
 * @param b 多项式系数
 * @param ab 结果
 * @param n 多项式长度
 * @param p 模数（质数）
 */
template <typename T>
void poly_multiply_ntt_cuda(T *a, T *b, T *ab, T n, T p, T OMEGA = 3);

/**
 * @brief 使用CUDA的NTT优化多项式乘法 - 优化版本v2（智能共享内存）
 *
 * 这个版本在v1基础上增加了智能共享内存优化，显著提升性能
 * 对小mid层使用共享内存，大mid层使用全局内存，避免数据依赖问题
 *
 * @param a 多项式系数
 * @param b 多项式系数
 * @param ab 结果
 * @param n 多项式长度
 * @param p 模数（质数）
 */
void poly_multiply_ntt_cuda_optimized_v2(u64 *a, u64 *b, u64 *ab, u64 n, u64 p,
                                         u64 OMEGA = 3);

/**
 * @brief 使用CUDA的NTT优化多项式乘法 - 优化版本v3（智能共享内存）
 *
 * 这个版本在v2基础上增加了智能共享内存优化，显著提升性能
 * 对小mid层使用共享内存，大mid层使用全局内存，避免数据依赖问题
 *
 * @param a 多项式系数
 * @param b 多项式系数
 * @param ab 结果
 * @param n 多项式长度
 * @param p 模数（质数）
 */
void poly_multiply_ntt_cuda_optimized_v3(u64 *a, u64 *b, u64 *ab, u64 n, u64 p,
                                         u64 OMEGA = 3);

/**
 * @brief GPU多项式乘法v1版本（预计算旋转因子优化）
 */
void poly_multiply_gpu_optimized_v1(u64 *a, u64 *b, u64 *result, u64 n, u64 p);

/**
 * @brief GPU多项式乘法v2版本（共享内存优化）
 */
void poly_multiply_gpu_optimized_v2(u64 *a, u64 *b, u64 *result, u64 n, u64 p);

/**
 * @brief GPU多项式乘法v2修复版本（安全共享内存优化）
 */
void poly_multiply_gpu_optimized_v2_fixed(u64 *a, u64 *b, u64 *result, u64 n, u64 p);

/**
 * @brief GPU多项式乘法v2增强版本（性能优化）
 */
void poly_multiply_gpu_optimized_v2_enhanced(u64 *a, u64 *b, u64 *result, u64 n, u64 p);

/**
 * @brief GPU多项式乘法v3版本（多层共享内存实验）
 */
void poly_multiply_gpu_optimized_v3_test3(u64 *a, u64 *b, u64 *result, u64 n, u64 p);

// CUDA kernel函数声明
template <typename T>
__global__ void ntt_forward_kernel(T *a_mont, T n, T p, T omega_mont, T mid);

template <typename T>
__global__ void ntt_inverse_kernel(T *a_mont, T n, T p, T omega_mont, T mid);

template <typename T>
__global__ void pointwise_mul_kernel(T *a_mont, T *b_mont, T *ab_mont, T n,
                                     T p);

// 辅助函数声明
template <typename T>
void check_cuda_error(const char *msg);

template <typename T>
void allocate_gpu_memory(T **ptr, size_t size);

template <typename T>
void copy_to_gpu(T *gpu_ptr, T *cpu_ptr, size_t size);

template <typename T>
void copy_to_cpu(T *cpu_ptr, T *gpu_ptr, size_t size);

template <typename T>
void free_gpu_memory(T *ptr);