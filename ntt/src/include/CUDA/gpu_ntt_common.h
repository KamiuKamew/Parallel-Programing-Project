#pragma once
#include "memory_safe_wrapper.h"
#include <iostream>

// 类型定义
typedef int64_t s64;

// 常量定义
const u64 NTT_P = 7340033;
const u64 NTT_OMEGA = 3;

// 基础模运算函数 - 内联定义
__device__ __host__ inline u64 mod_add(u64 a, u64 b, u64 p)
{
    u64 result = a + b;
    return result >= p ? result - p : result;
}

__device__ __host__ inline u64 mod_sub(u64 a, u64 b, u64 p)
{
    return a >= b ? a - b : a + p - b;
}

__device__ __host__ inline u64 mod_mul(u64 a, u64 b, u64 p)
{
    return ((unsigned __int128)a * b) % p;
}

__host__ inline u64 mod_inv(u64 a, u64 p)
{
    s64 x = 0, last_x = 1;
    s64 y = 1, last_y = 0;
    s64 r = p, last_r = a;

    while (r != 0)
    {
        s64 quotient = last_r / r;

        s64 temp = r;
        r = last_r - quotient * r;
        last_r = temp;

        temp = x;
        x = last_x - quotient * x;
        last_x = temp;

        temp = y;
        y = last_y - quotient * y;
        last_y = temp;
    }

    return last_x < 0 ? last_x + p : last_x;
}

__host__ inline u64 mod_pow(u64 base, u64 exp, u64 p)
{
    u64 result = 1;
    base %= p;
    while (exp > 0)
    {
        if (exp & 1)
            result = mod_mul(result, base, p);
        exp >>= 1;
        base = mod_mul(base, base, p);
    }
    return result;
}

inline void syncAndCheck(const char *operation)
{
    cudaError_t err = cudaDeviceSynchronize();
    if (err != cudaSuccess)
    {
        std::cout << "❌ CUDA错误 [" << operation
                  << "]: " << cudaGetErrorString(err) << std::endl;
        throw std::runtime_error("CUDA kernel执行失败");
    }
}