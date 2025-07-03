#include <iostream>
#include <vector>
#include <cstdint>

using u64 = uint64_t;

// 声明GPU NTT函数
void poly_multiply_gpu_optimized_v2(u64 *a, u64 *b, u64 *result, u64 n);

int main()
{
    std::cout << "=== GPU NTT简单测试 ===" << std::endl;

    // 测试数据：n=4的简单多项式乘法
    const u64 n = 4;
    std::vector<u64> a = {1, 2, 3, 4};
    std::vector<u64> b = {1, 2, 3, 4};
    std::vector<u64> result(2 * n - 1, 0);

    std::cout << "输入多项式 a: [";
    for (u64 i = 0; i < n; i++)
    {
        std::cout << a[i];
        if (i < n - 1)
            std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    std::cout << "输入多项式 b: [";
    for (u64 i = 0; i < n; i++)
    {
        std::cout << b[i];
        if (i < n - 1)
            std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    // 调用GPU NTT多项式乘法
    try
    {
        poly_multiply_gpu_optimized_v2(a.data(), b.data(), result.data(), n);

        std::cout << "GPU计算结果: [";
        for (u64 i = 0; i < 2 * n - 1; i++)
        {
            std::cout << result[i];
            if (i < 2 * n - 2)
                std::cout << ", ";
        }
        std::cout << "]" << std::endl;

        std::cout << "✅ GPU NTT测试完成！" << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cout << "❌ 测试失败: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}