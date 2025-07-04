#pragma once

#include <vector>
#include <string>

// 声明参考代码的主要接口函数
std::vector<long long> multiply_ntt_gpu(std::vector<long long> &poly1, std::vector<long long> &poly2,
                                        long long mod, long long primitive_root, const std::string &method = "naive");

// 包装函数，兼容原有接口
template <typename T>
inline void poly_multiply_ntt_gpu(T *a, T *b, T *ab, T n, T p, const std::string &method = "naive")
{
    using namespace std;

    // 将输入拷贝到 std::vector
    vector<long long> vec_a(n), vec_b(n);
    for (T i = 0; i < n; ++i)
    {
        vec_a[i] = static_cast<long long>(a[i]);
        vec_b[i] = static_cast<long long>(b[i]);
    }

    // 调用GPU实现
    std::vector<long long> result = multiply_ntt_gpu(vec_a, vec_b, static_cast<long long>(p), 3, method);

    // 将结果写回
    size_t target_len = result.size();
    for (size_t i = 0; i < target_len; ++i)
    {
        ab[i] = static_cast<T>(result[i]);
    }
}