/*
注：CRT部分仅支持T=u64, T2=u128，因此不进行模板化。
*/

#pragma once

#include "../general/utils.h"
#include "../OpenMP_Barrett/ntt.h"
#include "../CRT/const.h"
#include "../CRT/crt.h"
#include <algorithm>

#include <mpi.h>
#ifndef USE_MPI
#define USE_MPI
#endif

/**
 * @brief 使用NTT和CRT优化的多项式乘法 (MPI多进程版本)
 *
 * @param a 多项式系数 (assumed non-negative)
 * @param b 多项式系数 (assumed non-negative)
 * @param ab 结果多项式系数, modulo p
 * @param n 多项式长度 (number of coefficients, e.g., degree n-1)
 * @param p 最终结果的模数
 */
inline void poly_multiply_ntt_mpi(u64 *a, u64 *b, u64 *ab, u64 n, u64 p)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 广播基本参数
    MPI_Bcast((void *)&n, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast((void *)&p, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

    u64 n_expanded = expand_n(2 * n - 1);

    // 广播输入数据到所有进程
    MPI_Bcast((void *)a, n, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast((void *)b, n, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast((void *)&n_expanded, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

    // 只有主进程需要分配完整的CRT结果数组
    u64 **ab_crt = nullptr;
    u128 *ab_u128 = nullptr;

    if (rank == 0)
    {
        ab_crt = new u64 *[CRT_NUMS];
        for (u64 i = 0; i < CRT_NUMS; i++)
        {
            ab_crt[i] = new u64[n_expanded]{};
        }
        ab_u128 = new u128[n_expanded];
    }

    // 每个进程负责的CRT模数数量
    u64 mods_per_process = (CRT_NUMS + size - 1) / size; // 向上取整
    u64 start_mod = rank * mods_per_process;
    u64 end_mod = std::min(start_mod + mods_per_process, (u64)CRT_NUMS);

    // 为当前进程分配本地结果数组
    u64 local_mod_count = end_mod - start_mod;
    u64 **local_ab_crt = nullptr;
    if (local_mod_count > 0)
    {
        local_ab_crt = new u64 *[local_mod_count];
        for (u64 i = 0; i < local_mod_count; i++)
        {
            local_ab_crt[i] = new u64[n_expanded]{};
        }
    }

    // 每个进程计算分配给它的CRT模数
    for (u64 i = 0; i < local_mod_count; i++)
    {
        u64 mod_idx = start_mod + i;

        // 创建临时数组，对输入系数进行模数归约
        u64 *a_mod = new u64[n];
        u64 *b_mod = new u64[n];
        for (u64 j = 0; j < n; j++)
        {
            a_mod[j] = a[j] % CRT_MODS[mod_idx];
            b_mod[j] = b[j] % CRT_MODS[mod_idx];
        }

        poly_multiply_ntt_omp_Barrett(a_mod, b_mod, local_ab_crt[i], n,
                                      CRT_MODS[mod_idx], CRT_ROOTS[mod_idx]);

        delete[] a_mod;
        delete[] b_mod;
    }

    // 收集所有结果到主进程
    if (rank == 0)
    {
        // 主进程：复制自己的结果
        for (u64 i = 0; i < local_mod_count; i++)
        {
            for (u64 j = 0; j < n_expanded; j++)
            {
                ab_crt[start_mod + i][j] = local_ab_crt[i][j];
            }
        }

        // 接收其他进程的结果
        for (int src_rank = 1; src_rank < size; src_rank++)
        {
            u64 src_start = src_rank * mods_per_process;
            u64 src_end = std::min(src_start + mods_per_process, (u64)CRT_NUMS);
            u64 src_count = src_end - src_start;

            if (src_count > 0)
            {
                for (u64 i = 0; i < src_count; i++)
                {
                    MPI_Recv(ab_crt[src_start + i], n_expanded, MPI_UNSIGNED_LONG_LONG,
                             src_rank, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

        // 执行CRT合并
        for (u64 i = 0; i < n_expanded; ++i)
            ab_u128[i] = ab_crt[0][i];

        CRT_combine(ab_u128, ab_crt, n_expanded);

        // 最终模数归约
        for (u64 i = 0; i < n_expanded; ++i)
            ab[i] = ab_u128[i] % p;
    }
    else
    {
        // 非主进程：发送结果给主进程
        for (u64 i = 0; i < local_mod_count; i++)
        {
            MPI_Send(local_ab_crt[i], n_expanded, MPI_UNSIGNED_LONG_LONG,
                     0, i, MPI_COMM_WORLD);
        }
    }

    // 将最终结果广播给所有进程
    MPI_Bcast((void *)ab, n_expanded, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

    // 清理内存
    if (local_ab_crt)
    {
        for (u64 i = 0; i < local_mod_count; i++)
        {
            delete[] local_ab_crt[i];
        }
        delete[] local_ab_crt;
    }

    if (rank == 0)
    {
        if (ab_crt)
        {
            for (u64 i = 0; i < CRT_NUMS; i++)
            {
                delete[] ab_crt[i];
            }
            delete[] ab_crt;
        }
        delete[] ab_u128;
    }
}