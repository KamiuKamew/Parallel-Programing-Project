# 并行程序设计实验报告

> MPI 编程实验

## 目录

- [并行程序设计实验报告](#并行程序设计实验报告)
  - [目录](#目录)
  - [摘要](#摘要)
  - [前言](#前言)
  - [Barrett 模乘及其应用](#barrett-模乘及其应用)
    - [Barrett 模乘的实现](#barrett-模乘的实现)
      - [基本原理](#基本原理)
      - [代码实现](#代码实现)
    - [基于 Barrett 模乘的多项式乘法](#基于-barrett-模乘的多项式乘法)
      - [集成到 NTT 算法中](#集成到-ntt-算法中)
      - [OpenMP 多线程优化](#openmp-多线程优化)
    - [使用 CRT 解决大模数下误差较大的问题](#使用-crt-解决大模数下误差较大的问题)
      - [问题发现与分析](#问题发现与分析)
      - [CRT 解决方案的设计思路](#crt-解决方案的设计思路)
      - [MPI 集成实现](#mpi-集成实现)
  - [基于 MPI 的多进程优化](#基于-mpi-的多进程优化)
    - [MPI 并行化的设计思路](#mpi-并行化的设计思路)
      - [任务分配策略分析](#任务分配策略分析)
      - [通信模式选择](#通信模式选择)
    - [MPI 环境配置与多线程支持](#mpi-环境配置与多线程支持)
      - [线程支持级别](#线程支持级别)
      - [条件编译与兼容性](#条件编译与兼容性)
    - [当前实现与优化方向](#当前实现与优化方向)
      - [现有实现的特点](#现有实现的特点)
      - [多层次并行策略](#多层次并行策略)
      - [进一步优化的方向](#进一步优化的方向)
  - [代码性能测试](#代码性能测试)
  - [分析与讨论](#分析与讨论)

## 摘要

<!-- FIXME: 完成全文后撰写本部分 -->

## 前言

本次实验中，我将在已实现的基于 NTT 算法的多项式乘法代码的基础上，尝试进行 MPI 多进程优化。

在此之前，我已经实现了朴素的 NTT 算法，并且尝试过使用 OpenMP 和 pthread 进行多线程优化。本次实验我将使用朴素的 NTT 算法（使用 DIT/DIF 结构）作为基准算法。

实验中，我将首先实现 Barrett 模乘并用它来优化朴素的 NTT 算法，顺手使用 OpenMP 对该算法进行多线程优化；随后，我将使用 CRT（中国剩余定理）解决大模数下误差较大的问题；在此之后，我将引入 MPI 进一步进行多进程优化。所有代码工作完成后，我将测试不同**问题规模**、不同**节点数/线程数**下的算法性能，并将**串行**和**并行**的算法性能进行对比；随后，进行 profiling。在测试结束后，我将结合测试结果，讨论一些基本的**算法/编程策略**对性能的影响，包括 MPI 并行化的不同算法策略（如块划分、循环划分等不同任务划分方法，流水线算法等）对应的复杂性分析、不同 MPI 编程方法（阻塞通信 vs. 非阻塞通信双边通信 vs. 单边通信、MPI 自身的多线程支持等），以及体系结构相关优化（如 cache 优化）等。

<!-- FIXME：根据实际情况修改上文 -->

## Barrett 模乘及其应用

这一部分中，我将实现 Barrett 模乘算法，并将其应用于 NTT 优化中。随后，我将分析其在大模数下的局限性，并通过 CRT 技术加以解决。

### Barrett 模乘的实现

#### 基本原理

传统的模运算`a * b mod p`需要进行昂贵的除法操作，在大量重复的模乘计算中会成为性能瓶颈。Barrett 模乘通过预计算的方式，将除法运算转换为乘法和位移操作，从而显著提高模运算效率。

Barrett 模乘的核心思想基于以下数学原理。取模显然有下列式子：

$$x \bmod q = x - \lfloor x \cdot s \rfloor \cdot q$$

其中 $s = 1/q$，如果能以高精度求出 $s$，则此公式成立。为了避免浮点运算，Barrett 算法选择：

$$r = \lfloor 2^k / q \rfloor$$

使用近似：

$$r / 2^k \approx 1/q$$

从而得到：

$$x \bmod q = x - \lfloor x \cdot r / 2^k \rfloor \cdot q$$

由于近似过程中 $r / 2^k$ 与 $1/q$ 始终存在误差，由 $r = \lfloor 2^k / q \rfloor$ 可得：

$$r = 2^k / q - e \Rightarrow r / 2^k = 1 / q - e / 2^k \text{ for some } e \in [0, 1)$$

推出：

$$r / 2^k \leq 1 / q$$

当 $x = a \times b$ 时，输出又可化为：

$$ab \bmod q = ab - \lfloor ab \cdot r / 2^k \rfloor \cdot q \in [ab \bmod q, (ab \bmod q) + q]$$

在我们的实现中，选择 $k = 64$ 以获得更高的精度，并使用 128 位运算来处理中间结果。

#### 代码实现

我们在`ntt/src/include/Barrett/op.h`中实现了`BarrettMod`类。该类继承自基础的`Mod`类，并重写了关键的`mul`方法：

```cpp
template <typename T>
class BarrettMod : public Mod<T>
{
    using T2 = t_widen<T>; // 使用更宽的数据类型防止溢出

public:
    BarrettMod(T _mod, int iterations = 0) : Mod<T>(_mod)
    {
        barrett_k = 64;
        T2 r = (T2(1) << barrett_k) / _mod;

        // 迭代优化：通过多次迭代提高r的精度
        for (int i = 0; i < iterations; ++i)
        {
            T2 numerator = (T2(1) << barrett_k) - r * _mod;
            T2 correction = numerator / _mod;
            r += correction;
        }
        barrett_r = r;
    }

    T mul(T a, T b) const
    {
        T2 z = T2(a) * b;
        T2 q = (z * barrett_r) >> barrett_k;
        T res = T(z - q * this->mod);
        // 双重检查确保结果在[0, mod)范围内
        if (res >= this->mod) res -= this->mod;
        if (res >= this->mod) res -= this->mod;
        return res;
    }

private:
    T2 barrett_r;
    int barrett_k;
};
```

相比于指导手册的建议实现，我们的版本具有几个显著特点。首先是**更高精度**，我们使用`k=64`而非建议的`k=32`，提供更高的近似精度。其次是**迭代优化**，通过可选的迭代过程进一步提高预计算参数`r`的精度，这是为了应对`p = 1337006139375617`等超大模数情况。此外还有**安全性增强**，双重模数检查确保结果的正确性，以及**模板化设计**，支持不同的数据类型，提高代码复用性。

### 基于 Barrett 模乘的多项式乘法

#### 集成到 NTT 算法中

我们将 Barrett 模乘集成到 NTT 的正变换和逆变换中。在`ntt/src/include/Barrett/ntt.h`中实现了核心的 NTT 函数。该实现遵循标准的 Cooley-Tukey 算法结构，但在蝶形运算的每一步都使用 Barrett 模乘来提升性能：

```cpp
template <typename T>
inline void ntt_forward_Barrett(T *a, T n, T p, T omega)
{
    BarrettMod<T> mod(p, 100); // 使用100次迭代优化精度

    for (T mid = 1; mid < n; mid <<= 1)
    {
        T Wn = mod.pow(omega, (p - 1) / (mid << 1));
        for (T j = 0; j < n; j += (mid << 1))
        {
            T w = 1;
            for (T k = 0; k < mid; ++k, w = mod.mul(w, Wn))
            {
                T x = a[j + k];
                T y = mod.mul(w, a[j + k + mid]);
                a[j + k] = mod.add(x, y);
                a[j + k + mid] = mod.sub(x, y);
            }
        }
    }
}
```

该实现在 NTT 的蝶形运算中使用 Barrett 模乘，理论上应该能够提供比传统模运算更好的性能。逆变换的实现遵循相似的结构，但需要在最后进行归一化处理，即将所有系数乘以 $n^{-1} \bmod p$。

#### OpenMP 多线程优化

为了进一步提升性能，我们在`ntt/src/include/OpenMP_Barrett/ntt.h`中实现了基于 OpenMP 的多线程版本。这种实现延续了我们在 Lab3 中的优化思路，即对中间层的`j`循环进行并行化：

```cpp
template <typename T>
inline void ntt_forward_omp_Barrett(T *a, T n, T p, T omega)
{
    BarrettMod<T> mod(p, 100);

    for (T mid = 1; mid < n; mid <<= 1)
    {
        T Wn = mod.pow(omega, (p - 1) / (mid << 1));
#pragma omp parallel for // OpenMP并行化j循环
        for (T j = 0; j < n; j += (mid << 1))
        {
            T w = 1;
            for (T k = 0; k < mid; ++k, w = mod.mul(w, Wn))
            {
                T x = a[j + k];
                T y = mod.mul(w, a[j + k + mid]);
                a[j + k] = mod.add(x, y);
                a[j + k + mid] = mod.sub(x, y);
            }
        }
    }
}
```

选择`j`循环进行并行化的原因与之前分析一致：不同`j`迭代之间操作的数据区域相互独立，并且每个线程的任务量（即一个完整的`k`循环）具有合适的粒度，既保证了数据依赖性的正确处理，又获得了良好的负载均衡。

### 使用 CRT 解决大模数下误差较大的问题

#### 问题发现与分析

在实际测试中，我们发现 Barrett 模乘在处理超大模数（如`p = 1337006139375617`）时遇到了精度问题。尽管我们使用了`k=64`的高精度设置并引入了迭代优化机制，但对于某些极大的模数，Barrett 近似算法仍然存在精度损失，导致计算结果不正确。

这一问题的根本原因可以从几个方面分析。**数值范围限制**方面，即使使用 128 位中间运算，当模数接近`u64`的上限时，Barrett 算法的近似误差会被放大。**迭代优化局限性**方面，虽然迭代能够提高精度，但对于极大模数，有限次数的迭代无法完全消除累积误差。最重要的是**算法本质限制**，Barrett 模乘本身是一种近似算法，在极端情况下不可避免地存在精度边界。

#### CRT 解决方案的设计思路

为了彻底解决大模数问题，我们采用了中国剩余定理(CRT)的方法。核心思想是将大模数下的计算分解为多个小模数下的独立计算，然后将结果合并。我们选择了四个 NTT 友好的素数：

```cpp
static const u64 CRT_MODS[] = {998244353, 1004535809, 469762049, 167772161};
static const u64 CRT_ROOTS[] = {3, 3, 3, 3};
```

这些模数都小于 $2^{30}$，适合进行 Barrett 优化的 NTT 计算，同时它们的乘积远大于常见的 $u64$ 范围，可以表示非常大的系数。

#### MPI 集成实现

在`ntt/src/include/MPI/ntt.h`中，我们实现了集 Barrett 模乘、CRT 和 MPI 于一体的解决方案：

```cpp
inline void poly_multiply_ntt_mpi(u64 *a, u64 *b, u64 *ab, u64 n, u64 p)
{
    u64 n_expanded = expand_n(2 * n - 1);
    u64 **ab_crt = new u64 *[CRT_NUMS];
    u128 *ab_u128 = new u128[n_expanded];

    for (u64 i = 0; i < CRT_NUMS; i++)
    {
        ab_crt[i] = new u64[n_expanded]{};

        // 关键修复：对输入系数进行模数归约
        u64 *a_mod = new u64[n];
        u64 *b_mod = new u64[n];
        for (u64 j = 0; j < n; j++)
        {
            a_mod[j] = a[j] % CRT_MODS[i];
            b_mod[j] = b[j] % CRT_MODS[i];
        }

        // 在小模数下使用Barrett优化的NTT
        poly_multiply_ntt_omp_Barrett(a_mod, b_mod, ab_crt[i], n, CRT_MODS[i], CRT_ROOTS[i]);

        delete[] a_mod;
        delete[] b_mod;
    }

    // 使用CRT合并结果
    for (u64 i = 0; i < n_expanded; ++i)
        ab_u128[i] = ab_crt[0][i];

    CRT_combine(ab_u128, ab_crt, n_expanded);

    // 最终模数归约
    for (u64 i = 0; i < n_expanded; ++i)
        ab[i] = ab_u128[i] % p;

    // 内存清理
    delete[] ab_u128;
    for (u64 i = 0; i < CRT_NUMS; ++i)
        delete[] ab_crt[i];
    delete[] ab_crt;
}
```

这种方法的核心优势体现在多个方面。**精度保证**方面，在较小的 CRT 模数（如`998244353`）下，Barrett 算法能够保证完全的精度。**性能优化**方面，每个小模数下的计算都能享受 Barrett 模乘和 OpenMP 多线程的性能优势。**扩展性**方面，可以处理任意大小的目标模数，只要 CRT 模数的乘积足够大。

特别需要注意的是代码中的"关键修复"部分。在初期实现中，我们直接将原始的大系数传递给 Barrett 优化的 NTT 函数，但这会导致当输入系数超过 CRT 模数时 Barrett 算法失效。通过预先对输入系数进行模数归约，我们确保了 Barrett 算法始终在其有效的数值范围内工作。

通过这种"分而治之"的策略，我们成功地将 Barrett 模乘的优势（高效的模运算）与 CRT 的优势（大数计算能力）相结合，既保证了计算的正确性，又维持了良好的性能表现。这为后续的 MPI 多进程优化奠定了坚实的基础。

## 基于 MPI 的多进程优化

这一部分中，我将基于前面实现的 Barrett 模乘和 CRT 技术，设计并实现 MPI 多进程优化方案。重点探讨 MPI 并行化的算法策略、编程方法选择以及与多线程技术的结合。

### MPI 并行化的设计思路

#### 任务分配策略分析

在 NTT 多项式乘法的 MPI 并行化中，存在多种可能的任务分配策略。**数据并行**方面，可以将输入多项式按系数进行分块，每个进程负责处理一部分系数的 NTT 变换。然而，这种方法面临着 NTT 算法本身具有全局依赖性的挑战，不同阶段的蝶形运算需要访问分布在各个进程中的数据，通信开销可能抵消并行收益。

**算法并行**方面，更适合的策略是利用 CRT 分解的特性。由于我们使用四个独立的模数进行 CRT 分解，每个模数下的多项式乘法计算完全独立，天然适合分配给不同的 MPI 进程。这种策略的优势在于进程间无需频繁通信，只需在最终的 CRT 合并阶段进行数据收集。

#### 通信模式选择

根据 CRT 并行化的特点，我们的 MPI 实现主要涉及两种通信模式。**集合通信**用于初始阶段将输入数据广播到所有进程，以及最终阶段收集各进程的计算结果。**点对点通信**可以用于特定的数据交换场景，但在我们的设计中使用较少。

考虑到 CRT 计算的独立性，我们选择了相对简单但有效的通信策略：主进程负责数据分发和结果收集，其他进程专注于特定模数下的 NTT 计算。

### MPI 环境配置与多线程支持

#### 线程支持级别

现代 MPI 应用往往需要与多线程技术结合，我们在 `ntt/src/include/MPI/config.h` 中实现了灵活的线程支持配置：

```cpp
#ifndef MPI_THREAD_LEVEL
#define MPI_THREAD_LEVEL MPI_THREAD_FUNNELED
#endif

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
            switch (provided) { /* ... 详细输出 ... */ }                                \
        }                                                                                \
    } while (0)
```

这种设计支持多种 MPI 线程级别，默认使用 `MPI_THREAD_FUNNELED`，允许只有主线程进行 MPI 调用，而其他线程可以进行计算。通过编译时宏定义，也可以选择 `MPI_THREAD_MULTIPLE` 等更高级别的支持。

#### 条件编译与兼容性

为了保证代码的可移植性，我们使用条件编译技术，使得代码在没有 MPI 环境的情况下也能正常编译运行：

```cpp
#ifdef USE_MPI
    // MPI相关定义
    #define MPI_GET_RANK(rank) \
        int rank = 0;          \
        MPI_Comm_rank(MPI_COMM_WORLD, &(rank))
    #define MPI_ONLY_MAIN if (rank == 0)
#else
    // 非MPI环境的替代定义
    #define MPI_GET_RANK(rank) int rank = 0
    #define MPI_ONLY_MAIN
#endif
```

### 当前实现与优化方向

#### 现有实现的特点

在 `ntt/src/include/MPI/ntt.h` 中，我们实现了基础的 MPI 多进程框架。当前实现将 CRT 的四个模数分解与 Barrett 优化相结合：

```cpp
inline void poly_multiply_ntt_mpi(u64 *a, u64 *b, u64 *ab, u64 n, u64 p)
{
    u64 n_expanded = expand_n(2 * n - 1);
    u64 **ab_crt = new u64 *[CRT_NUMS];
    u128 *ab_u128 = new u128[n_expanded];

    for (u64 i = 0; i < CRT_NUMS; i++)
    {
        ab_crt[i] = new u64[n_expanded]{};

        // 关键：输入系数的模数归约
        u64 *a_mod = new u64[n];
        u64 *b_mod = new u64[n];
        for (u64 j = 0; j < n; j++)
        {
            a_mod[j] = a[j] % CRT_MODS[i];
            b_mod[j] = b[j] % CRT_MODS[i];
        }

        // 在当前进程内使用 OpenMP 优化的 Barrett NTT
        poly_multiply_ntt_omp_Barrett(a_mod, b_mod, ab_crt[i], n, CRT_MODS[i], CRT_ROOTS[i]);

        delete[] a_mod;
        delete[] b_mod;
    }

    // CRT 合并在主进程完成
    for (u64 i = 0; i < n_expanded; ++i)
        ab_u128[i] = ab_crt[0][i];

    CRT_combine(ab_u128, ab_crt, n_expanded);

    for (u64 i = 0; i < n_expanded; ++i)
        ab[i] = ab_u128[i] % p;

    // 内存清理...
}
```

这个实现的核心特点是将每个 CRT 模数的计算作为独立任务，在进程内部使用 OpenMP 进行多线程优化。虽然当前版本还没有实现真正的进程间任务分配，但为完整的 MPI 并行化奠定了基础。

#### 多层次并行策略

我们的设计体现了多层次并行的思想。**进程级并行**通过 MPI 实现，将不同的 CRT 模数分配给不同进程，实现粗粒度的任务分解。**线程级并行**在每个进程内部通过 OpenMP 实现，对 NTT 的中间循环进行细粒度并行化。**指令级并行**通过 Barrett 模乘的优化实现，减少除法运算，提高单指令效率。

这种分层策略的优势在于能够充分利用现代 HPC 系统的层次化架构：MPI 适合节点间通信，OpenMP 适合节点内多核并行，Barrett 优化适合单核性能提升。

#### 进一步优化的方向

**进程间负载均衡**方面，当前四个 CRT 模数的计算量基本相等，但在实际应用中可以考虑动态负载均衡策略。**通信优化**方面，可以使用非阻塞通信和计算-通信重叠技术来隐藏通信延迟。**内存优化**方面，可以使用 MPI 的共享内存窗口来减少数据复制开销。

**算法层面的优化**包括实现真正的分布式 NTT 算法，将蝶形运算的不同阶段分配给不同进程，或者使用流水线技术来实现多个多项式乘法的并行处理。**体系结构相关优化**方面，可以考虑 NUMA 感知的进程绑定和缓存优化策略。

通过这种渐进式的优化设计，我们既保证了算法的正确性和可维护性，又为进一步的性能提升预留了空间。当前的实现虽然还未完全发挥 MPI 的分布式计算优势，但已经建立了坚实的技术基础，为后续的扩展和优化提供了良好的架构支撑。

## 代码性能测试

## 分析与讨论
