# 并行程序设计实验报告

> 基于 pthread 和 OpenMP 的 NTT 算法多线程优化

## 目录

- [并行程序设计实验报告](#并行程序设计实验报告)
  - [目录](#目录)
  - [前言](#前言)
  - [基于 OpenMP 的朴素算法多线程优化](#基于-openmp-的朴素算法多线程优化)
    - [对优化位点的理论分析](#对优化位点的理论分析)
    - [进行基于 OpenMP 的多线程并行化](#进行基于-openmp-的多线程并行化)
    - [初步的性能测试与分析](#初步的性能测试与分析)
  - [基于 pthread 的朴素算法多线程优化](#基于-pthread-的朴素算法多线程优化)
  - [CRT 优化算法的实现](#crt-优化算法的实现)
    - [1. 基本原理](#1-基本原理)
    - [2. 代码实现](#2-代码实现)
  - [基于 pthread 的 CRT 优化算法多线程优化](#基于-pthread-的-crt-优化算法多线程优化)
  - [分析与讨论](#分析与讨论)

## 前言

本次实验中，我将在已实现的基于 NTT 算法的多项式乘法代码的基础上，尝试进行多线程优化。

在此之前，我已经在上一次实验中实现了朴素的 NTT 算法以及使用 Montgomery 规约优化的 NTT 算法。本次实验我将选择使用 Montgomery 规约优化的 NTT 算法作为基准算法。

实验中，我将首先使用 OpenMP 和 pthread 在基准算法上进行多线程优化；随后，我将使用 CRT（中国剩余定理）对基准算法进行优化，并在优化后的算法的基础上使用 pthread 进行多线程优化。最后，我将测试不同**问题规模**、不同**线程数**下的算法性能（串行和并行对比）；讨论一些基本的**算法/编程策略**对性能的影响，以及 **Pthread 程序和 OpenMP 程序**的**性能差异**；讨论多线程并行化的**不同算法策略**（如矩阵水平划分、垂直划分等不同任务划分方法，不同算法策略下的一致性保证等、线程管理代价优化等）及其**复杂性分析**；**profiling** 及**体系结构相关优化**（如 cache 优化）；不同**平台**（x86 或 ARM）上并行化实验；**OpenMP 卸载到加速器设备**；与 oneAPI 编程、C++语言标准中的多线程编程进行**性能对比**等。

<!-- FIXME：这里提到的不全是要做的。最后根据做了什么修改一下即可 -->

本次实验的主要优化对象是`ntt_forward_mont`和`ntt_inverse_mont`两个函数。由于`ntt_forward_mont`和`ntt_inverse_mont`均由结构相似的三层`for`循环构成，本文将以`ntt_forward_mont`为例进行分析，`ntt_inverse_mont`的部分同理可知，因此略去不做展示。

## 基于 OpenMP 的朴素算法多线程优化

这一部分中，我将使用 OpenMP 在基准算法上进行多线程优化。

### 对优化位点的理论分析

`ntt_forward_mont` 函数的核心循环结构如下所示：

```cpp
// 第一层循环，记为 mid 循环
for (u32 mid = 1; mid < n; mid <<= 1)
{
    u32_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1));
    // 第二层循环，记为 j 循环
    for (u32 j = 0; j < n; j += (mid << 1))
    {
        u32_mont w_mont = montMod.from_u32(1); // 为每个 j 块初始化 w_mont
        // 第三层循环，记为 k 循环
        for (u32 k = 0; k < mid; ++k, w_mont = montMod.mul(w_mont, Wn_mont))
        {
            u32_mont x_mont = a_mont[j + k];
            u32_mont y_mont = montMod.mul(w_mont, a_mont[j + k + mid]);
            a_mont[j + k] = montMod.add(x_mont, y_mont);
            a_mont[j + k + mid] = montMod.sub(x_mont, y_mont);
        }
    }
}
```

为了确认并行化位点，我们首先对这三层循环的多线程并行化潜力进行分析：

1. 外层 `mid` 循环：

   - 此循环控制 NTT 算法的计算阶段。每个阶段（`mid` 的一次迭代）都**依赖于**前一阶段对数组 `a_mont` 的完整更新结果。由于 `a_mont` 是**就地修改**的，当前 `mid` 值的计算**必须在**前一个 `mid` 值的计算**完成后**进行。
   - 由于这种迭代间的强数据依赖，`mid` 循环本身是**不适合直接多线程并行化**的，其执行必须是串行的。

2. 中间 `j` 循环：

   - 对于固定的 `mid` 值，`j` 循环以 `mid << 1` 为步长遍历数据块。在 `j` 循环的每次迭代中，所操作的 `a_mont` 数组的索引范围（例如，`j1` 对应的 `[j1, j1 + (mid << 1) - 1]` 和 `j2` 对应的 `[j2, j2 + (mid << 1) - 1]`）是**完全不重叠**的。旋转因子 `w_mont` 在每个 `j` 迭代开始时被分别初始化为 `montMod.from_u32(1)`。
   - 由于不同 `j` 迭代之间操作的数据区域**相互独立**，并且迭代内的状态（如 `w_mont`）被分别初始化，因此 `j` 循环的各个迭代是**数据独立**的，**非常适合多线程并行化**。
   - 若并行化 `j` 循环，每个线程将负责执行其内部的整个 `k` 循环，即进行 `mid` 次蝶形运算。应注意到，线程的任务量（也就是任务粒度）会随着 `mid` 的增加（即算法阶段的深入）而相应增大（并不会出现粒度过细的情况）。不过，对于确定的`mid`，各个线程的**任务量几乎相同**，它们会大致同时完成。
   - 注：此处还可以进一步进行若干讨论，例如线程数量、反复创建线程带来的性能损耗，等等。不过这里只是实现一个朴素版本，因此此处不做进一步讨论，而是留到[基于 OpenMP 的朴素算法多线程优化](#基于-openmp-的朴素算法多线程优化)这一节中进行统一分析。

3. 内层 `k` 循环

   - 在 `k` 循环内部，对 `a_mont[j + k]` 和 `a_mont[j + k + mid]` 的写操作，对于不同的 `k` 值，访问的是不同的内存单元，因此**不存在写冲突**。然而，旋转因子 `w_mont` 是通过 `w_mont = montMod.mul(w_mont, Wn_mont)` **迭代更新**的。这意味着第 `k` 次迭代所使用的 `w_mont` 值依赖于第 `k-1` 次迭代更新后的 `w_mont` 值，这构成了**循环携带依赖**。
   - 如果**简单地**在 `k` 循环前（即 L2 循环内部）应用并行化指令，会遇到 `w_mont` 的**依赖问题**。同时，在每个 `j` 的迭代中都创建和销毁线程组将导致巨大的并行开销。**不过**理论上，可以通过修改 `w_mont` 的计算方式来**解除此依赖**，例如在每个 `k` 迭代中独立计算 `w_mont_for_k = montMod.pow(Wn_mont, k)`（考虑到初始 `w_mont` 的值）。若**如此修改**，`k` 循环的迭代**便可独立执行**。
   - 即便在解决了 `w_mont` 的依赖后并行化 `k` 循环，每个并行任务也仅执行一次蝶形运算。这属于**非常细的任务粒度**。这会导致任务过度划分，降低整体性能。

综上，针对《2024 并行程序设计 Lab3_1_Pthread 编程-1.pdf》中提到的"显然可以对**第三层循环**进行多线程优化"，虽然理论上第三层循环在处理好 `w_mont` 的迭代依赖（例如通过为每个 `k` 重新计算 `(W_n)^k`）后可以并行，但这种做法存在任务粒度过细、并行开销高、实现难度大的问题。相比之下，并行化第二层循环（即 `j` 循环）任务粒度合适、并行管理开销较低、实现简洁。

因此，综合考虑数据依赖性、任务粒度、并行开销和实现复杂度，我们选择对**第二层** `for` 循环（`j` 循环）进行并行化。

### 进行基于 OpenMP 的多线程并行化

为了使用 OpenMP 对第二层 `for` 循环进行并行化，我们只需要在 `ntt_forward_mont` 函数中对应的 `j` 循环前添加 OpenMP 的并行化指导语句 `#pragma omp parallel for`（如前所述，`ntt_inverse_mont`同理）。这样，编译器会在运行时自动将循环迭代分配给多个线程执行。即：

```cpp
// ...
    for (u32 mid = 1; mid < n; mid <<= 1)
    {
        u32_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1));
#pragma omp parallel for // 添加此行
        for (u32 j = 0; j < n; j += (mid << 1))
        {
            u32_mont w_mont = montMod.from_u32(1);
            for (u32 k = 0; k < mid; ++k, w_mont = montMod.mul(w_mont, Wn_mont))
            {
                // ... 运算主体 ...
            }
        }
    }
// ...
```

### 初步的性能测试与分析

为了编译启用了 OpenMP 的代码，我们使用`g++ -std=c++17 -O2 -fopenmp -I<include_path> -o ntt_parallel ntt_test_program.cpp`进行编译，并以类似的指令编译了一个串行版本`ntt_serial`（不带`-fopenmp`）。我们设计了一个测试脚本，该脚本针对多个预设的测试用例（标记为 0.in, 1.in, 2.in, 3.in，其问题规模各有不同）运行这些程序。为了减少性能波动带来的误差，每个测试配置（特定测试用例+特定线程数）均重复执行 5 次，并记录平均、最小和最大执行时间。测试用例`4.in`由于其模数超出了程序当前`u32`类型的处理范围，在本次初步测试中被跳过。

初步测试验证了 OpenMP 优化后的代码的**正确性**，各测试样例均能正确通过。针对小规模问题(即`0.in`)和中大规模问题(即`1.in`，`2.in`和`3.in`)，测试结果有较大区别：

- 在小规模问题上，OpenMP 并行版本（即使是单线程）的平均执行时间远高于纯串行版本。串行执行时间约为 3 **微秒**，而单线程 OpenMP 平均耗时接近 2 **毫秒**，16 线程版本甚至达到平均 9 毫秒。这符合预期，即对于任务粒度极小的计算，并行化引入的额外开销（线程创建、管理、同步等）会超过并行执行本身带来的潜在收益。
- 在中大规模问题上，OpenMP 表现出了一定的**加速**效果。通常从 2 个线程开始，平均执行时间相较于串行版本有所降低。例如，在`2.in`上，串行平均 111 毫秒，2 线程 OpenMP 平均约 99 毫秒，8 线程时平均约 71 毫秒。最佳加速效果通常出现在 4 至 12 个线程之间，具体取决于测试用例。然而，加速比并非随线程数线性增长，且在线程数较多（如 16 线程）时，部分测试用例性能**略有下降**。

我们将在报告的后续"分析与讨论"章节中，进一步设计测试方法，得到更全面的性能数据，从而进行更深入和系统的分析。

<!-- FIXME：根据后文实际撰写内容修改这段话内容 -->

## 基于 pthread 的朴素算法多线程优化

这一部分中，我将使用 pthread 在基准算法上进行多线程优化。

我们的优化在上一节（[基于 OpenMP 的朴素算法多线程优化](#基于-openmp-的朴素算法多线程优化)）的基础上进行。经过分析，我们选择了中间 `j` 循环作为优化位点。然而，由于该层循环并非最外层循环，如果我们在中间 `j` 循环外（即外层 `mid` 循环内）创建线程，则会导致线程被反复创建与删除，大大降低加速效果。因此，在这一节中，我们将抛弃 OpenMP，使用 pthread 维护一个线程池，动态地向中间 `j` 循环分配线程。

<!-- TODO -->

## CRT 优化算法的实现

这一部分中，我将使用 CRT（中国剩余定理）对基准算法进行优化。

### 1. 基本原理

在使用 NTT 进行多项式乘法时，若乘积多项式的系数可能非常大，以至于超过单个 NTT 模数 `p` 所能表示的范围（或者为了提高 NTT 的效率，选用的 `p` 较小），直接使用单一模数进行计算将导致结果错误。中国剩余定理（CRT）为此提供了一种解决方案，其核心思想是将一个大数上的计算分解为在多个较小的、互质的模数上进行计算，然后将这些结果合并以获得原始大数域上的解。

<!-- 或许这里可以加一下CRT的数学原理。 -->

在多项式乘法的背景下，这意味着：

1. 选择一组素数 `m_0, m_1, ..., m_{k-1}`，这些素数都适合进行 NTT（即 `m_i - 1` 具有足够大的 2 的幂次因子），并且它们的乘积 `M = m_0 * m_1 * ... * m_{k-1}` 必须大于多项式乘积结果的任何可能系数的最大值。
2. 对于每个模数 `m_i`，独立地计算多项式乘积 `C_i(x) = A(x) * B(x) (mod m_i)`。这通常通过对输入多项式 `A(x)` 和 `B(x)` 的系数分别取模 `m_i`，然后执行标准的 NTT 乘法完成。
3. 对于结果多项式的每一个系数，我们得到一组同余方程：
   `c_j ≡ c_{j,0} (mod m_0)`
   `c_j ≡ c_{j,1} (mod m_1)`
   `...`
   `c_j ≡ c_{j,k-1} (mod m_{k-1})`
   其中 `c_j` 是最终结果多项式第 `j` 个系数，`c_{j,i}` 是在模 `m_i` 下计算得到的第 `j` 个系数。
4. 使用 CRT 从 `c_{j,0}, c_{j,1}, ..., c_{j,k-1}` 解出 `c_j (mod M)`。

### 2. 代码实现

在项目 `ntt/src/include/CRT/ntt.h` 文件中，我们实现了基于 CRT 的 NTT 多项式乘法。

首先，定义了一组预选的 NTT 友好素数 `CRT_MODS` 及其对应的原根 `CRT_ROOTS`：

```cpp
static const u64 CRT_MODS[] = {998244353, 1004535809, 469762049, 167772161};
static const u64 CRT_ROOTS[] = {3, 3, 3, 3};
static const u64 CRT_NUMS = sizeof(CRT_MODS) / sizeof(CRT_MODS[0]);
```

这些模数都小于 `2^30`，适合 `u64` 计算，并且它们的乘积远大于常见的 `u64` 范围，可以表示非常大的系数。

核心函数 `poly_multiply_ntt_crt` 负责整个流程：

```cpp
inline void poly_multiply_ntt_crt(u64 *a, u64 *b, u64 *ab, u64 n, u64 p)
{
    u64 n_expanded = expand_n(2 * n - 1); // 确定NTT运算的长度

    u64 **ab_crt = new u64 *[CRT_NUMS];    // 存储每个模数下的NTT结果
    u128 *ab_u128 = new u128[n_expanded]; // 存储CRT合并后的结果 (使用u128防止溢出)

    // 1. 对每个CRT模数执行标准NTT多项式乘法
    for (u64 i = 0; i < CRT_NUMS; i++)
    {
        ab_crt[i] = new u64[n_expanded]{};
        // poly_multiply_ntt 对 a 和 b 的系数模 CRT_MODS[i] 后进行NTT乘法
        poly_multiply_ntt(a, b, ab_crt[i], n, CRT_MODS[i], CRT_ROOTS[i]);
    }

    // 初始化合并结果，以第一个模数下的结果为基础
    for (u64 i = 0; i < n_expanded; ++i)
        ab_u128[i] = ab_crt[0][i];

    // 2. 使用CRT合并结果
    CRT_combine(ab_u128, ab_crt, n_expanded);
    // 或者 CRT_combine_2(ab_u128, ab_crt, n_expanded);

    // 3. 将合并后的大数结果对最终模数 p 取模
    for (u64 i = 0; i < n_expanded; ++i)
        ab[i] = ab_u128[i] % p;

    // 内存回收
    delete[] ab_u128;
    for (u64 i = 0; i < CRT_NUMS; ++i)
        delete[] ab_crt[i];
    delete[] ab_crt;
}
```

该函数首先对选定的每个 `CRT_MODS[i]` 执行一次标准的多项式乘法 `poly_multiply_ntt`，并将结果存储在 `ab_crt[i]` 中。随后，调用 `CRT_combine` (或 `CRT_combine_2`) 函数，依据中国剩余定理，将这些部分结果合并到 `ab_u128` 数组中（使用 `u128` 类型以容纳可能的大数值）。最后，将合并后的结果对目标模数 `p` 取模，得到最终的多项式系数。

代码中实现了两种 CRT 合并算法：`CRT_combine` 和 `CRT_combine_2`。
`CRT_combine` 的实现方式如下（一种迭代式的 Garner\'s algorithm 的变体）：

```cpp
// 将 ab_crt 的 CRT 结果合并到 ab 中 (ab 初始化为 ab_crt[0])
inline void CRT_combine(u128 *ab, u64 **ab_crt, u64 n)
{
    u128 m = CRT_MODS[0]; // 当前已经合并的模数的乘积
    for (u64 i = 1; i < CRT_NUMS; ++i) // 从第二个模数开始迭代
    {
        u64 CRT_MOD = CRT_MODS[i]; // 当前要合并的模数
        Mod128 mod(CRT_MOD);       // 用于在 CRT_MOD 下进行运算的模运算类

        // 计算 m 在模 CRT_MOD 下的逆元：inv(m) mod CRT_MOD
        u128 inv = mod.inv(m % CRT_MOD);
        for (u64 j = 0; j < n; j++) // 对每个系数进行合并
        {
            u128 x = ab[j]; // 当前已经合并的结果 x = k_0 (mod m_0*...*m_{i-1})
            // t = (ab_crt[i][j] - (x % CRT_MOD)) * inv(m) (mod CRT_MOD)
            // t 是满足 (x + m*t) === ab_crt[i][j] (mod CRT_MOD) 的最小非负整数
            u64 t = mod.sub(ab_crt[i][j], x % CRT_MOD);
            t = mod.mul(t, inv);

            // 更新合并结果：x_new = x + m * t
            // 此时 x_new 满足 x_new === k_0 (mod m_0*...*m_{i-1})
            // 且 x_new === ab_crt[i][j] (mod CRT_MOD)
            x = x + m * t;
            ab[j] = x;
        }
        m *= CRT_MOD; // 更新已合并模数的乘积
    }
}
```

此方法逐个引入新的模数，并更新当前已合并的解。对于每个系数 `ab[j]`，它首先保存了模 `m_0 * ... * m_{i-1}` 的解，然后通过计算一个调整项 `t`，使得新的解也满足模 `m_i` 的同余条件。

`CRT_combine_2` 实现了另一种合并方式，它对每个系数独立地从头开始构建 CRT 解：

```cpp
inline void CRT_combine_2(u128 *ab, u64 **ab_crt, u64 n)
{
    for (u64 i = 0; i < n; ++i) // 对每个系数独立计算
    {
        u128 x = ab_crt[0][i]; // 从第一个模数的结果开始
        u128 m = CRT_MODS[0];

        for (u64 j = 1; j < CRT_NUMS; ++j) // 逐个引入其他模数的结果
        {
            u64 CRT_MOD = CRT_MODS[j];
            Mod128 mod(CRT_MOD);

            // t = (ab_crt[j][i] - (x % CRT_MOD)) * inv(m % CRT_MOD) (mod CRT_MOD)
            u64 t = mod.sub(ab_crt[j][i], x % CRT_MOD);
            u64 inv = mod.inv(m % CRT_MOD);
            t = mod.mul(t, inv);

            x = x + m * t; // 更新解
            m *= CRT_MOD;  // 更新模数乘积
        }
        ab[i] = x; // 存储当前系数的最终CRT解
    }
}
```

我们分别对两个 CRT 算法在第 0 到第 4 个测试样例上进行初步的正确性测试和性能测试。初步结果显示，两个算法均可以正确解决这 5 个测试样例，证明算法在`uint64`的数据宽度下也可以正常工作； `CRT_combine` 在处理 `n = 131072` 规模的问题时，耗时约 **440-450** 微秒，而 `CRT_combine_2` 耗时约 **900-910** 微秒。这表明迭代更新当前解的 `CRT_combine` 方法可能具有**更好的缓存局部性**或**更少的重复计算**（例如模逆元的计算 `mod.inv(m % CRT_MOD)` 中，`m` 在 `CRT_combine` 的外层循环中，而 `CRT_combine_2` 中 `m` 在内层循环变化，但 `mod.inv` 每次都针对新的 `m % CRT_MOD` 进行计算）。我们将在[分析与讨论](#分析与讨论)一节中对该现象的成因进行分析。

引入这一 CRT 优化算法后，我们可以处理更大系数的多项式乘法。我们将在接下来的[基于 pthread 的 CRT 优化算法多线程优化](#基于-pthread-的-crt-优化算法多线程优化)一节中对这一算法进行并行化。

<!--
注：使用第一个CRT算法的结果如下：

多项式乘法结果正确
average latency for n = 4 p = 7340033 : 0.03446 (us)
多项式乘法结果正确
average latency for n = 131072 p = 7340033 : 446.686 (us)
多项式乘法结果正确
average latency for n = 131072 p = 104857601 : 442.015 (us)
多项式乘法结果正确
average latency for n = 131072 p = 469762049 : 441.801 (us)
多项式乘法结果正确
average latency for n = 131072 p = 1337006139375617 : 448.786 (us)

使用第二个CRT算法的结果如下：

多项式乘法结果正确
average latency for n = 4 p = 7340033 : 0.04449 (us)
多项式乘法结果正确
average latency for n = 131072 p = 7340033 : 905.087 (us)
多项式乘法结果正确
average latency for n = 131072 p = 104857601 : 905.88 (us)
多项式乘法结果正确
average latency for n = 131072 p = 469762049 : 900.773 (us)
多项式乘法结果正确
average latency for n = 131072 p = 1337006139375617 : 910.255 (us)

我们可以在这里简单提一下，并在分析与讨论中的"2. "或者"6. "下进行讨论。
 -->

## 基于 pthread 的 CRT 优化算法多线程优化

这一部分中，我将在优化后的算法的基础上使用 pthread 进行多线程优化。

由于我们的模数数量固定为 4 个，我们可以简单地把四次多项式乘法分配给四个线程。此时线程划分次数少，每个线程的任务量几乎一致，且创建线程次数少，因此没必要使用线程池，直接为四个多项式乘法分配四个线程即可。

我们主要修改了 `poly_multiply_ntt_crt` 函数（在 pthread 实现中更名为 `poly_multiply_ntt_pthread_crt`），并引入了辅助的数据结构和线程工作函数。其核心思路是将针对 `CRT_NUMS`（固定为 4）个不同模数的 `poly_multiply_ntt` 调用分配给不同的线程并行执行。

我们定义了一个结构体 `PthreadNttArgs`，用于封装传递给每个线程的独立参数。

```cpp
// Structure to pass arguments to each NTT worker thread
struct PthreadNttArgs
{
    u64 *a_poly;          // Pointer to the first input polynomial
    u64 *b_poly;          // Pointer to the second input polynomial
    u64 *result_poly_crt; // Pointer to the output array for this thread (a part of ab_crt)
    u64 n_poly_len;       // Original length of the polynomials
    u64 current_mod;      // The CRT modulus for this thread
    u64 current_root;     // The primitive root for the current_mod
};
```

每个被创建的线程将执行 `poly_multiply_ntt_thread_worker` 函数。此函数从传入的参数中解析出所需数据，并调用标准的 `poly_multiply_ntt` 函数完成特定模数下的多项式乘法。

```cpp
// Thread worker function: performs poly_multiply_ntt for a single CRT modulus
static void *poly_multiply_ntt_thread_worker(void *arg)
{
    PthreadNttArgs *params = (PthreadNttArgs *)arg;
    poly_multiply_ntt(params->a_poly, params->b_poly, params->result_poly_crt,
                        params->n_poly_len, params->current_mod, params->current_root);
    return NULL;
}
```

在主函数 `poly_multiply_ntt_pthread_crt` 中，我们首先为每个模数的结果数组分配内存。然后，创建 `CRT_NUMS` 个线程，每个线程配置其独立的 `PthreadNttArgs`。通过 `pthread_create` 启动这些线程后，主线程通过 `pthread_join` 等待所有子线程完成计算。

```cpp
inline void poly_multiply_ntt_pthread_crt(u64 *a, u64 *b, u64 *ab, u64 n, u64 p)
{
    u64 n_expanded = expand_n(2 * n - 1);

    u64 **ab_crt = new u64 *[CRT_NUMS];
    u128 *ab_u128 = new u128[n_expanded];

    // Step 1: Allocate memory for each CRT result array (serially)
    for (u64 i = 0; i < CRT_NUMS; i++)
    {
        ab_crt[i] = new u64[n_expanded]{};
    }

    pthread_t threads[CRT_NUMS];
    PthreadNttArgs thread_args[CRT_NUMS];

    // Step 2: Create and launch threads to perform NTT for each modulus in parallel
    for (u64 i = 0; i < CRT_NUMS; i++)
    {
        thread_args[i].a_poly = a;
        thread_args[i].b_poly = b;
        thread_args[i].result_poly_crt = ab_crt[i];
        thread_args[i].n_poly_len = n;
        thread_args[i].current_mod = CRT_MODS[i];
        thread_args[i].current_root = CRT_ROOTS[i];

        pthread_create(&threads[i], NULL, poly_multiply_ntt_thread_worker, (void *)&thread_args[i]);
    }

    // Step 3: Wait for all threads to complete their execution
    for (u64 i = 0; i < CRT_NUMS; i++)
    {
        pthread_join(threads[i], NULL);
    }

    // Initialize ab_u128 with results from the first modulus
    for (u64 i = 0; i < n_expanded; ++i)
        ab_u128[i] = ab_crt[0][i];

    // Combine results using CRT
    CRT_combine(ab_u128, ab_crt, n_expanded);

    // Final modular reduction
    for (u64 i = 0; i < n_expanded; ++i)
        ab[i] = ab_u128[i] % p;

    // Memory cleanup
    delete[] ab_u128;
    for (u64 i = 0; i < CRT_NUMS; ++i)
        delete[] ab_crt[i];
    delete[] ab_crt;
}
```

我们在`./test.sh 2 1 1`和`./test.sh 2 4 4`下分别 pthread 优化后的 CRT 进行初步的正确性测试和性能测试。初步结果显示，优化后的代码可以正确解决这 5 个测试样例，证明算法在 pthread 优化后也可以正常工作；二者对规模为`n = 131072`的样例的处理时长均从未优化时的$440 \mu s$变为$140 \mu s$左右，一方面证明 pthread 优化确实有效，另一方面似乎表明`test.sh`里设置的核心数和线程数似乎与程序实际使用的线程数无关。

<!--
测试结果：

在./test.sh 2 1 1（2表示测试pthread，1表示申请1个核心，1表示申请1个线程）下，结果如下：

多项式乘法结果正确
average latency for n = 4 p = 7340033 : 0.327471 (us)
多项式乘法结果正确
average latency for n = 131072 p = 7340033 : 140.988 (us)
多项式乘法结果正确
average latency for n = 131072 p = 104857601 : 141.119 (us)
多项式乘法结果正确
average latency for n = 131072 p = 469762049 : 141.79 (us)
多项式乘法结果正确
average latency for n = 131072 p = 1337006139375617 : 145.227 (us)

然而，在./test.sh 2 4 4下仍然差不多是以上结果。
 -->

## 分析与讨论

这一部分中，我将尝试进行以下分析：

1. 测试不同**问题规模**、不同**线程数**下的算法性能（串行和并行对比）
2. 讨论一些基本的**算法/编程策略**对性能的影响
3. 讨论 **Pthread 程序和 OpenMP 程序**的**性能差异**
4. 讨论多线程并行化的**不同算法策略**，及其**复杂性分析**
   1. 例如：矩阵水平划分、垂直划分等不同任务划分方法，不同算法策略下的一致性保证等、线程管理代价优化等
5. 进行**profiling**
6. 进行**体系结构相关优化**（如 cache 优化）；
7. 在不同**平台**（x86 或 ARM）上进行并行化实验
8. 尝试**OpenMP 卸载到加速器设备**
9. 与 oneAPI 编程、C++语言标准中的多线程编程进行**性能对比**。

<!-- FIXME：这里提到的不全是要做的。最后根据做了什么修改一下即可 -->
