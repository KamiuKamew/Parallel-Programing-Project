#pragma once

#ifdef USE_MPI

#ifndef MPI_THREAD_LEVEL
#define MPI_THREAD_LEVEL MPI_THREAD_FUNNELED
#endif
// 也可以使用 mpic++ -O3 -fopenmp -DUSE_MPI -DMPI_THREAD_LEVEL=MPI_THREAD_MULTIPLE -o main main.cc 编译来指定线程级别

#include <mpi.h>
// #define MPI_INIT(argc, argv)        \
//     do                              \
//     {                               \
//         MPI_Init(&(argc), &(argv)); \
//     } while (0)
// #define MPI_INIT(argc, argv)                                            \
//     do                                                                  \
//     {                                                                   \
//         int provided;                                                   \
//         MPI_Init_thread(&(argc), &(argv), MPI_THREAD_LEVEL, &provided); \
//     } while (0)
#define MPI_INIT(argc, argv)                                                             \
    do                                                                                   \
    {                                                                                    \
        int provided;                                                                    \
        MPI_Init_thread(&(argc), &(argv), MPI_THREAD_LEVEL, &provided);                  \
        if (provided < MPI_THREAD_LEVEL)                                                 \
        {                                                                                \
            std::cerr << "警告：MPI 提供的线程级别 (" << provided                        \
                      << ") 低于请求的等级 (" << MPI_THREAD_LEVEL << ")！" << std::endl; \
        }                                                                                \
        else                                                                             \
        {                                                                                \
            std::cout << "[MPI] 线程支持等级初始化成功：";                               \
            switch (provided)                                                            \
            {                                                                            \
            case MPI_THREAD_SINGLE:                                                      \
                std::cout << "MPI_THREAD_SINGLE\n";                                      \
                break;                                                                   \
            case MPI_THREAD_FUNNELED:                                                    \
                std::cout << "MPI_THREAD_FUNNELED\n";                                    \
                break;                                                                   \
            case MPI_THREAD_SERIALIZED:                                                  \
                std::cout << "MPI_THREAD_SERIALIZED\n";                                  \
                break;                                                                   \
            case MPI_THREAD_MULTIPLE:                                                    \
                std::cout << "MPI_THREAD_MULTIPLE\n";                                    \
                break;                                                                   \
            default:                                                                     \
                std::cout << "未知等级：" << provided << "\n";                           \
                break;                                                                   \
            }                                                                            \
        }                                                                                \
    } while (0)
#define MPI_FINALIZE()  \
    do                  \
    {                   \
        MPI_Finalize(); \
    } while (0)
#define MPI_GET_RANK(rank) \
    int rank = 0;          \
    MPI_Comm_rank(MPI_COMM_WORLD, &(rank))
#define MPI_ONLY_MAIN if (rank == 0)
#define MPI_BARRIER() MPI_Barrier(MPI_COMM_WORLD)

// MPI 计时宏定义
#define TIMER_START() double start_time = MPI_Wtime()
#define TIMER_END() double end_time = MPI_Wtime()
#define TIMER_ELAPSED() ((end_time - start_time) * 1000.0) // 转换为毫秒

#else

#define MPI_INIT(argc, argv)
#define MPI_FINALIZE()
#define MPI_GET_RANK(rank) int rank = 0
#define MPI_ONLY_MAIN
#define MPI_BARRIER()

// 标准 chrono 计时宏定义
#include <chrono>
#define TIMER_START() auto start_time = std::chrono::high_resolution_clock::now()
#define TIMER_END() auto end_time = std::chrono::high_resolution_clock::now()
#define TIMER_ELAPSED() (std::chrono::duration<double, std::ratio<1, 1000>>(end_time - start_time).count()) // 毫秒

#endif

#ifdef USE_CUDA

#define CUDA_WARMUP()      \
    int *dummy;            \
    cudaMalloc(&dummy, 1); \
    cudaFree(dummy);       \
    cudaDeviceSynchronize();

#else

#define CUDA_WARMUP()

#endif