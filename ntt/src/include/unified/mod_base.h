#pragma once

#include "../general/type.h"

/**
 * @brief 统一模运算基类
 *
 * 定义了所有模运算算法的通用接口，支持数域转换
 */
class ModBase
{
public:
    ModBase(u32 _mod) : mod(_mod) {}
    virtual ~ModBase() = default;

    // 禁用拷贝构造和赋值
    ModBase(const ModBase &) = delete;
    ModBase &operator=(const ModBase &) = delete;

    // 基本模运算接口
    virtual u32 add(u32 a, u32 b) const = 0;
    virtual u32 sub(u32 a, u32 b) const = 0;
    virtual u32 mul(u32 a, u32 b) const = 0;
    virtual u32 pow(u32 base, u32 exp) const = 0;
    virtual u32 inv(u32 x) const = 0;

    // 数域转换接口
    virtual u32 to_compute_domain(u32 a) const = 0;
    virtual u32 from_compute_domain(u32 a) const = 0;

    // 数组批量转换接口
    virtual void array_to_compute_domain(u32 *a, u32 n) const;
    virtual void array_from_compute_domain(u32 *a, u32 n) const;

    // 获取模数
    u32 get_mod() const { return mod; }

    // 获取算法名称（用于调试）
    virtual const char *get_algorithm_name() const = 0;

protected:
    u32 mod;
};

// 数组批量转换的默认实现
inline void ModBase::array_to_compute_domain(u32 *a, u32 n) const
{
    for (u32 i = 0; i < n; ++i)
        a[i] = to_compute_domain(a[i]);
}

inline void ModBase::array_from_compute_domain(u32 *a, u32 n) const
{
    for (u32 i = 0; i < n; ++i)
        a[i] = from_compute_domain(a[i]);
}