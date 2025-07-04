#pragma once

#include <vector>
#include <string>
#include <stdint.h>

#ifdef __NVCC__
#include <cuda_runtime.h>
#endif

using u32 = uint32_t;
using u64 = uint64_t;

// 基础32位模运算
__device__ __host__ inline u32 mul_mod(u32 a, u32 b, u32 mod)
{
    return (u64)a * b % mod;
}

__device__ __host__ inline u32 add_mod(u32 a, u32 b, u32 mod)
{
    u32 sum = a + b;
    return (sum >= mod) ? (sum - mod) : sum;
}

__device__ __host__ inline u32 sub_mod(u32 a, u32 b, u32 mod)
{
    return (a >= b) ? (a - b) : (a + mod - b);
}

__device__ __host__ inline u32 pow_mod(u32 base, u32 exp, u32 mod)
{
    u32 result = 1;
    base %= mod;
    while (exp > 0)
    {
        if (exp & 1)
            result = mul_mod(result, base, mod);
        base = mul_mod(base, base, mod);
        exp >>= 1;
    }
    return result;
}

__device__ __host__ inline u32 inv_mod(u32 a, u32 mod)
{
    return pow_mod(a, mod - 2, mod);
}

// Barrett Reducer for 32-bit
struct BarrettReducer32
{
    u32 mod;
    u64 mu;

    __host__ __device__ BarrettReducer32(u32 m) : mod(m)
    {
        if (m != 0)
        {
            mu = ((1ULL << 32) / m) + 1;
        }
        else
        {
            mu = 0;
        }
    }

    __host__ __device__ u32 multiply(u32 a, u32 b) const
    {
        if (mod == 0)
            return 0;
        u64 prod = (u64)a * b;
        u64 q = (prod * mu) >> 32;
        u32 r = (u32)(prod - q * mod);
        if (r >= mod)
            r -= mod;
        return r;
    }
};

// Montgomery Reducer for 32-bit
struct MontgomeryReducer32
{
    u32 mod;
    u32 neg_mod_inv;
    u32 r2;
    static constexpr int R_BITS = 32;
    static constexpr u64 R = 1ULL << R_BITS;

    __host__ __device__ MontgomeryReducer32(u32 m = 0) : mod(m), neg_mod_inv(0), r2(0)
    {
        if (mod == 0)
            return;

        // 计算 -mod^(-1) mod R
        u32 inv = 1;
        for (int i = 0; i < 5; ++i)
        { // 5次牛顿迭代足够32位
            inv = inv * (2 - mod * inv);
        }
        neg_mod_inv = -inv;

        // 计算 R^2 mod mod
        u64 r_mod = R % mod;
        r2 = (u32)((r_mod * r_mod) % mod);
    }

    __host__ __device__ u32 to_mont(u32 a) const
    {
        return multiply(a, r2);
    }

    __host__ __device__ u32 from_mont(u32 a_mont) const
    {
        u64 t = a_mont;
        u32 m = (u32)t * neg_mod_inv;
        u32 u = (u32)((t + (u64)m * mod) >> R_BITS);
        if (u >= mod)
            u -= mod;
        return u;
    }

    __host__ __device__ u32 multiply(u32 a_mont, u32 b_mont) const
    {
        u64 t = (u64)a_mont * b_mont;
        u32 m = (u32)t * neg_mod_inv;
        u32 u = (u32)((t + (u64)m * mod) >> R_BITS);
        if (u >= mod)
            u -= mod;
        return u;
    }
};

// GPU NTT函数声明
std::vector<u32> multiply_ntt_gpu32(std::vector<u32> &poly1, std::vector<u32> &poly2,
                                    u32 mod, u32 primitive_root, const std::string &method = "basic");

// 调试用：简单的CPU版本NTT（用于对比）
std::vector<u32> multiply_ntt_cpu32_debug(std::vector<u32> &poly1, std::vector<u32> &poly2,
                                          u32 mod, u32 primitive_root);

// 包装函数，兼容原有接口
template <typename T>
inline void poly_multiply_ntt_gpu32(T *a, T *b, T *ab, T n, T p, const std::string &method = "basic")
{
    using namespace std;

    // 将输入拷贝到 std::vector
    vector<u32> vec_a(n), vec_b(n);
    for (T i = 0; i < n; ++i)
    {
        vec_a[i] = static_cast<u32>(a[i]);
        vec_b[i] = static_cast<u32>(b[i]);
    }

    // 调用GPU实现
    std::vector<u32> result = multiply_ntt_gpu32(vec_a, vec_b, static_cast<u32>(p), 3, method);

    // 将结果写回
    size_t target_len = result.size();
    for (size_t i = 0; i < target_len; ++i)
    {
        ab[i] = static_cast<T>(result[i]);
    }
}