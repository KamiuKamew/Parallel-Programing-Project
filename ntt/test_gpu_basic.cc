#include "src/include/CUDA/ntt.h"
#include "src/include/transform.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <cstring>

// 测试参数
const u64 TEST_N = 8;       // 从小规模开始测试
const u64 TEST_P = 7340033; // 模数
const u64 TEST_OMEGA = 3;   // 原根

void print_array(const char *name, u64 *arr, u64 n)
{
    printf("%s: ", name);
    for (u64 i = 0; i < n; i++)
    {
        printf("%llu ", arr[i]);
    }
    printf("\n");
}

bool test_gpu_vs_cpu()
{
    printf("=== GPU vs CPU NTT测试 (n=%llu) ===\n", TEST_N);

    // 准备测试数据
    u64 a_cpu[TEST_N] = {1, 2, 3, 4, 0, 0, 0, 0};
    u64 a_gpu[TEST_N];
    memcpy(a_gpu, a_cpu, sizeof(a_cpu));

    printf("原始数据: ");
    print_array("", a_cpu, TEST_N);

    // CPU版本正变换
    printf("\n--- CPU版本 ---\n");
    ntt_forward(a_cpu, TEST_N, TEST_P, TEST_OMEGA);
    print_array("CPU NTT", a_cpu, TEST_N);

    // GPU版本正变换
    printf("\n--- GPU版本 ---\n");
    ntt_forward_gpu(a_gpu, TEST_N, TEST_P, TEST_OMEGA);
    print_array("GPU NTT", a_gpu, TEST_N);

    // 比较结果
    bool forward_match = true;
    for (u64 i = 0; i < TEST_N; i++)
    {
        if (a_cpu[i] != a_gpu[i])
        {
            printf("正变换不匹配: 位置%llu, CPU=%llu, GPU=%llu\n", i, a_cpu[i], a_gpu[i]);
            forward_match = false;
        }
    }

    if (forward_match)
    {
        printf("✓ 正变换结果匹配\n");
    }
    else
    {
        printf("✗ 正变换结果不匹配\n");
        return false;
    }

    return true;
}

int main()
{
    printf("开始GPU NTT基础功能测试...\n\n");

    // 检查CUDA设备
    int device_count;
    cudaError_t err = cudaGetDeviceCount(&device_count);
    if (err != cudaSuccess || device_count == 0)
    {
        printf("错误: 未检测到CUDA设备\n");
        return 1;
    }

    printf("检测到 %d 个CUDA设备\n\n", device_count);

    // 运行测试
    bool all_tests_passed = true;

    if (!test_gpu_vs_cpu())
    {
        all_tests_passed = false;
    }

    printf("\n=== 测试总结 ===\n");
    if (all_tests_passed)
    {
        printf("✓ 所有测试通过！GPU NTT实现正确。\n");
        return 0;
    }
    else
    {
        printf("✗ 部分测试失败，需要检查GPU实现。\n");
        return 1;
    }
}