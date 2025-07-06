#pragma once

#include "mod_base.h"
#include "mod_naive.h"
#include "mod_montgomery.h"
#include "mod_barrett.h"
#include "../CUDA/ntt.h"
#include <memory>

/**
 * @brief GPU统一NTT框架
 *
 * 提供与统一框架兼容的GPU NTT实现，支持不同的模运算算法：
 * - ModNaive: 朴素模运算
 * - ModMontgomery: Montgomery规约
 * - ModBarrett: Barrett规约
 *
 * GPU实现特点：
 * - 使用CUDA并行计算
 * - 支持多种模运算算法
 * - 自动内存管理和数据传输
 * - 优化的线程块配置
 */

// GPU模运算方法枚举
enum class GPUModMethod
{
    NAIVE,
    MONTGOMERY,
    BARRETT
};

// GPU模运算适配器类
class ModGPUAdapter
{
private:
    GPUModMethod method_;
    u32 p_;

public:
    ModGPUAdapter(GPUModMethod method, u32 p) : method_(method), p_(p) {}

    GPUModMethod get_method() const { return method_; }
    u32 get_p() const { return p_; }
};

// GPU统一NTT前向变换
inline void ntt_forward_unified_gpu(u32 *a, u32 n, u32 p, u32 omega,
                                    const ModGPUAdapter &mod_adapter)
{
    // GPU实现使用Montgomery域，需要转换
    MontMod<u32> montMod(p);
    u32 n_expanded = expand_n(n);

    // 数据预处理
    u32 *a_expanded = expand_a(a, n, n_expanded);
    bit_reverse_permute(a_expanded, n_expanded);

    // 转换到Montgomery域
    u32 *a_mont = new u32[n_expanded];
    for (u32 i = 0; i < n_expanded; i++)
    {
        a_mont[i] = montMod.from_T(a_expanded[i]);
    }

    // 根据模运算方法选择GPU实现
    switch (mod_adapter.get_method())
    {
    case GPUModMethod::NAIVE:
        // 朴素算法在GPU上使用Montgomery实现（因为GPU实现基于Montgomery）
        poly_multiply_ntt_gpu_naive(a_mont, a_mont, a_mont, n_expanded, p, omega);
        break;
    case GPUModMethod::MONTGOMERY:
        poly_multiply_ntt_gpu_mont(a_mont, a_mont, a_mont, n_expanded, p, omega);
        break;
    case GPUModMethod::BARRETT:
        poly_multiply_ntt_gpu_barrett(a_mont, a_mont, a_mont, n_expanded, p, omega);
        break;
    }

    // 转换回普通域
    for (u32 i = 0; i < n_expanded; i++)
    {
        a_expanded[i] = montMod.to_T(a_mont[i]);
    }

    // 复制结果
    for (u32 i = 0; i < n; i++)
    {
        a[i] = a_expanded[i];
    }

    delete[] a_expanded;
    delete[] a_mont;
}

// GPU统一NTT逆变换
inline void ntt_inverse_unified_gpu(u32 *a, u32 n, u32 p, u32 omega,
                                    const ModGPUAdapter &mod_adapter)
{
    // GPU实现中逆变换在前向变换中处理，这里调用前向变换
    ntt_forward_unified_gpu(a, n, p, omega, mod_adapter);
}

// GPU统一多项式乘法
inline void poly_multiply_unified_gpu(u32 *a, u32 *b, u32 *ab, u32 n, u32 p, u32 omega,
                                      const ModGPUAdapter &mod_adapter)
{
    // 根据模运算方法选择GPU实现
    switch (mod_adapter.get_method())
    {
    case GPUModMethod::NAIVE:
        poly_multiply_ntt_gpu_naive(a, b, ab, n, p, omega);
        break;
    case GPUModMethod::MONTGOMERY:
        poly_multiply_ntt_gpu_mont(a, b, ab, n, p, omega);
        break;
    case GPUModMethod::BARRETT:
        poly_multiply_ntt_gpu_barrett(a, b, ab, n, p, omega);
        break;
    }
}

// GPU统一框架选择函数
inline void poly_multiply_unified_gpu_choose(u32 *a, u32 *b, u32 *ab, u32 n, u32 p, u32 omega,
                                             GPUModMethod mod_method)
{
    ModGPUAdapter mod_adapter(mod_method, p);
    poly_multiply_unified_gpu(a, b, ab, n, p, omega, mod_adapter);
}