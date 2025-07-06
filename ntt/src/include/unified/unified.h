#pragma once

#include "mod_base.h"
#include "mod_naive.h"
#include "mod_montgomery.h"
#include "mod_barrett.h"
#include "ntt_unified.h"
#include "ntt_unified_openmp.h"
#include "ntt_unified_pthread.h"
// #include "ntt_unified_mpi.h"  // 暂时禁用MPI，需要特殊编译环境
#include "ntt_unified_gpu.h"
#include <memory>

/**
 * @brief 统一模运算框架
 *
 * 包含所有模运算算法的实现：
 * - ModNaive: 朴素模运算
 * - ModMontgomery: Montgomery规约
 * - ModBarrett: Barrett规约
 *
 * 以及统一的NTT算法实现：
 * - 串行版本: ntt_forward_unified, ntt_inverse_unified, poly_multiply_unified
 * - OpenMP版本: ntt_forward_unified_omp, ntt_inverse_unified_omp, poly_multiply_unified_omp
 * - Pthread版本: ntt_forward_unified_pthread, ntt_inverse_unified_pthread, poly_multiply_unified_pthread
 *
 * 所有算法都继承自ModBase，提供统一的接口
 */

enum class ModMethod
{
    NAIVE,
    MONTGOMERY,
    BARRETT
};

enum class ParallelMethod
{
    SERIAL,
    PTHREAD,
    OPENMP,
    MPI,
    GPU
};

inline void poly_multiply_unified_choosewith(u32 *a, u32 *b, u32 *ab, u32 n, u32 p, u32 omega, ModMethod mod_method, ParallelMethod parallel_method)
{
    std::shared_ptr<ModBase> mod_impl; // 为了兼容pthread的线程池实现，这里使用shared_ptr
    switch (mod_method)
    {
    case ModMethod::NAIVE:
        mod_impl = std::make_shared<ModNaive>(p);
        break;
    case ModMethod::MONTGOMERY:
        mod_impl = std::make_shared<ModMontgomery>(p);
        break;
    case ModMethod::BARRETT:
        mod_impl = std::make_shared<ModBarrett>(p);
        break;
    }

    switch (parallel_method)
    {
    case ParallelMethod::SERIAL:
        poly_multiply_unified(a, b, ab, n, p, omega, mod_impl.get());
        break;
    case ParallelMethod::PTHREAD:
        poly_multiply_unified_pthread(a, b, ab, n, p, omega, mod_impl.get());
        break;
    case ParallelMethod::OPENMP:
        poly_multiply_unified_omp(a, b, ab, n, p, omega, mod_impl.get());
        break;
    case ParallelMethod::MPI:
        poly_multiply_unified_mpi(a, b, ab, n, p, omega, mod_impl.get());
        break;
    case ParallelMethod::GPU:
        // 根据模运算方法选择GPU实现
        switch (mod_method)
        {
        case ModMethod::NAIVE:
            poly_multiply_ntt_gpu_naive(a, b, ab, n, p, omega);
            break;
        case ModMethod::MONTGOMERY:
            poly_multiply_ntt_gpu_mont(a, b, ab, n, p, omega);
            break;
        case ModMethod::BARRETT:
            poly_multiply_ntt_gpu_barrett(a, b, ab, n, p, omega);
            break;
        }
        break;
    }
}