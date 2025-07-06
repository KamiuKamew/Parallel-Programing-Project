#include "src/include/ntt.h"
#include "src/include/OpenMP/ntt.h"
#include "src/include/CRT/ntt.h"
#include "src/include/pthread_crt/ntt.h"
#include "src/include/pthread_simple/ntt.h"
#include "src/include/Barrett/ntt.h"
#include "src/include/OpenMP_Barrett/ntt.h"
#include "src/include/base4/ntt-base4.h"
#include "src/include/unified/ntt_unified_gpu.h"
#include "src/include/general/utils.h"
#include "src/include/config.h"

#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sys/time.h>
#include <chrono>

// GPU模式下的MPI处理
#ifdef USE_CUDA
#ifndef USE_MPI
// 在GPU模式下，如果没有启用MPI，则提供空的MPI宏
#define MPI_INIT(argc, argv)
#define MPI_FINALIZE()
#define MPI_GET_RANK(rank) int rank = 0
#define MPI_ONLY_MAIN
#define MPI_BARRIER()
#endif
#endif

std::string path_read = "/nttdata/";
std::string path_check = "/nttdata/";
std::string path_write = "files/";

template <typename T>
void fRead(T *a, T *b, T *n, T *p, T input_id)
{
    // 数据输入函数
    std::string str2 = std::to_string(input_id);
    std::string strin = path_read + str2 + ".in";
    char data_path[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), data_path);
    data_path[strin.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    fin >> *n >> *p;
    for (T i = 0; i < *n; i++)
    {
        fin >> a[i];
    }
    for (T i = 0; i < *n; i++)
    {
        fin >> b[i];
    }
}

template <typename T>
void fCheck(T *ab, T n, T input_id)
{
    // 判断多项式乘法结果是否正确
    std::string str2 = std::to_string(input_id);
    std::string strout = path_check + str2 + ".out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    for (T i = 0; i < n * 2 - 1; i++)
    {
        T x;
        fin >> x;
        if (x != ab[i])
        {
            std::cout << "多项式乘法结果错误" << std::endl;
            return;
        }
    }
    std::cout << "多项式乘法结果正确" << std::endl;
    return;
}

template <typename T>
void fWrite(T *ab, T n, T input_id)
{
    // 数据输出函数, 可以用来输出最终结果, 也可用于调试时输出中间数组
    std::string str2 = std::to_string(input_id);
    std::string strout = path_write + str2 + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (T i = 0; i < n * 2 - 1; i++)
    {
        fout << ab[i] << '\n';
    }
}

template <typename T>
void poly_multiply(T *a, T *b, T *ab, T n, T p)
{
    for (T i = 0; i < n; ++i)
    {
        for (T j = 0; j < n; ++j)
        {
            ab[i + j] = (1LL * a[i] * b[j] % p + ab[i + j]) % p;
        }
    }
}

template <typename T>
int _main(int argc, char *argv[])
{
    T a[300000], b[300000], ab[300000];

    MPI_GET_RANK(rank);

    // CUDA 预热，消除首次调用的冷启动开销
    CUDA_WARMUP();

    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 263882790666241
    // 第四个模数超过了整型表示范围, 如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求, 如果要自行探索大模数 NTT,
    // 请在完成前三个模数的基础代码及优化后实现大模数 NTT 输入文件共五个,
    // 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久,
    // 推荐调试正确性时只使用输入文件 1
    T test_begin = 0;
    T test_end = 3; // Skip test 4 which needs u64
    for (T i = test_begin; i <= test_end; ++i)
    {
        long double ans = 0;
        T n_, p_;

        MPI_ONLY_MAIN { fRead(a, b, &n_, &p_, i); }

        memset(ab, 0, sizeof(ab));
        TIMER_START();

        // GPU统一框架测试
        poly_multiply_unified_gpu_choose(a, b, ab, n_, p_, 3, GPUModMethod::MONTGOMERY);

        TIMER_END();
        ans += TIMER_ELAPSED();

        MPI_ONLY_MAIN
        {
            fCheck(ab, n_, i);
            std::cout << "average latency for n = " << n_ << " p = " << p_ << " : "
                      << ans << " (us) " << std::endl;
            // 可以使用 fWrite 函数将 ab 的输出结果打印到 files 文件夹下
            // 禁止使用 cout 一次性输出大量文件内容
            fWrite(ab, n_, i);
        }
    }
    return 0;
}

int main(int argc, char *argv[])
{

    MPI_INIT(argc, argv);

    path_read = "/home/hexay/projects/Lab Proj/Parallel-Programing-Project/.nttdata/";
    path_check = "/home/hexay/projects/Lab Proj/Parallel-Programing-Project/.nttdata/";
    path_write = "/home/hexay/projects/Lab Proj/Parallel-Programing-Project/.file/";

    int result = _main<u32>(argc, argv);

    MPI_FINALIZE();

    return result;
}