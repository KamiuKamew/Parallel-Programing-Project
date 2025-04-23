# 8 NTT 选题：SIMD 实验

NTT 的 SIMD 优化需要分别实现以下两个难点 (也即给分点):

1. neon 的**向量化操作**不支持**取模**，需要自行实现

2. NTT(递归版) 的循环主体需要使用蝴蝶变换，导致向量内部需要进行**手动的交换**

前提级，如果你实在无法实现某一个难点的向量化，也可以回退到最基础的版本 (朴素多项式乘法和逐个取模)，最终得分**根据**实现的 SIMD 优化**类型给分**，而非根据最终评测性能给分。

> 也就是说，只需要实现就行。

## 8.1 向量化取模的替代实现方案

**如果无法完成**向量化取模，甚至多项式乘法都难以实现，在这里提出三个解决方案:

### 8.1.1 取模操作串行化

涉及到取模时单独处理向量里的每一个元素，也就是除取模操作外其他操作均使用向量化。

此处对该解决方法不做描述。

### 8.1.2 浮点数近似取模

使用浮点数除法近似取模操作。浮点数近似取模原理比较简单:

$a \times b \% p = a \times b - p \times \left\lfloor \frac{a \times b}{p} \right\rfloor = a \times b - p \times \left\lfloor \frac{1}{p} \times a \times b \right\rfloor$

这一串式子中，$\frac{1}{p}$ 可以提前预处理出来，因此所有操作均可由 SIMD 提供的函数实现。

需要注意的是，由于模数较大，且涉及到浮点数操作，使用这种方法进行向量化取模会导致精度损失，进而导致答案与正确答案存在一定偏差，如果采用了这种方式实现 SIMD 取模**一定**会导致**答案错误**,但在本次实验中对这部分**误差**进行了**忽略**，即使你使用浮点数近似取模最终导致多项式乘法的答案错误(除样例 0 和样例 1 外甚至相差非常多)，也会**按照答案正确给分** (但不会超过 Montgomery 规约的分数)。

### 8.1.3 Montgomery 规约

使用 Montgomery 规约将模乘转化为支持向量化的操作。

在了解 Montgomery 规约前需要保证你对乘法逆元等基础数论概念有所了解 (即 NTT 前置数论知识)。

Montgomery 规约是专门用于取模优化的算法 (后续实验也会介绍另一种规约)，目前网络上介绍原理的资料较多，仅对原理进行简述，本次**实验核心**也**不在于** Montgomery 规约的**原理**而**是**如何**使用 SIMD 优化** Montgomery 规约，如果你**参考了**已有的某些网站 Montgomery 规约的代码或仓库，**必须在报告中指出**，否则会进行扣分 (当然 neon 优化的不可能找得到)，提醒，这一部分的实现难度相对较大。

Montgomery 空间由模数 n 和一个满足 $r \ge n$ 且与 n 互质的正整数 r 定义，通常取 $r = 2^{32}$ 或$r = 2^{64}$。

定义 $x$ 在 Montgomery 空间中的数值为 $\bar{x} = x \cdot r \mod n$。

在 Montgomery 空间内，加减等的运算与常规运算相同，但乘法不同，定义 $*$ 为 Montgomery 空间内的乘法，$\cdot$ 为常规运算乘法，则

$$
\bar{x} * \bar{y} = x \cdot y = (x \cdot y) \cdot r \mod n
$$

$$
\bar{x} \cdot \bar{y} = (x \cdot y) \cdot r \cdot r \mod n
$$

$$
\bar{x} * \bar{y} = \bar{x} \cdot \bar{y} \cdot r^{-1} \mod n
$$

对于 $\bar{x} \cdot \bar{y} \mod n$ 可以直接计算，对于涉及到除法 $x \cdot r^{-1} \mod n$ ，Montgomery 空间内除 2 不需要除法，因此可以得到最终化简式 ($n'$ 可求):

$$
x \cdot r^{-1} \equiv \frac{x - x \cdot n' \mod r \cdot n}{r}
$$

通过规约将乘法取模转化。  
一个简要的规约实现如下:

```cpp
typedef __uint32_t u32;
typedef __uint64_t u64;

const u32 n = 1e9 + 7，nr = inverse(n，1ull << 32);

u32 reduce(u64 x) {
    u32 q = u32(x) * nr; // q = x * n’ mod r
    u64 m = (u64) q * n; // m = q * n
    u32 y = (x - m) >> 32; // y = (x - m) / r
    return x < m ? y + n : y; // 保证 y 非负
}
```

容易发现规约中的所有操作均可以使用 SIMD 优化，你可以所有运算过程均在 Montgomery 数域下进行，也可以只针对模乘。

- [维基百科对 Montgomery 模乘的详细介绍](https://en.wikipedia.org/wiki/Montgomery_modular_multiplication)
- [上古时代介绍 Montgomery 规约原理的论文](https://www.ams.org/journals/mcom/1985-44-170/S0025-5718-1985-0777282-X/S0025-5718-1985-0777282-X.pdf)
- [知乎一篇介绍几种优化取模方法](https://zhuanlan.zhihu.com/p/16156667678)

## 8.2 向量化蝴蝶变换

### 8.2.1 基本实现框架

```cpp
for(int mid = 1; mid < limit; mid <<= 1) {
    for(int j = 0; j < limit; j += (mid << 1)) {
        int w = 1;// 旋转因子
        for( int k = 0; k < mid; k++, w = w * Wn) {
            // 运算主体
            // 计算 a[j + k],
            // 计算 a[j + k + mid];
        }
    }
}
```

以上是蝴蝶变换的主体，mid 代表当前运算的步长，第一层循环是对分治的模拟，第二层循环和第三层循环是具体的分治过程，即模拟合并两块等步长的序列。第三层循环显然可以进行 SIMD 优化，当步长小于向量长度时，直接使用朴素方法计算，当步长大于等于向量长度，就可以使用向量进行优化。

### 8.2.2 进阶实现：DIT/DIF

在完成以上蝴蝶变换的基础要求后，在这里提出一个较难的算法：DIT/DIF，其参考资料放在最下方。

- DIT(Decimation in Time)，按时间抽取的 NTT 的一种实现就是 Cooley-Tukey 算法，即常见的实现，这里直接跳过。
- DIF(Decimation in Frequency)，按频率抽取的 NTT。

DIT 的作用是输入一个正常顺序序列的 $a$，一个按位翻转后的序列 $\omega_{rev}$，输出按位翻转后的按时间抽取的信号。

DIF 的作用是输入一个按位翻转后的信号，一个按位翻转后的序列 $\omega^{-1}_{rev}$，输出正常顺序序列 $a$。所以只需要将多项式先 DIT，做完运算后直接 DIF，这样就不需要位翻转，只需要在初始化的时候进行一次。如果你尝试了这种方式优化，需要注意向量内部的变换。

- [DIF/DIT 的实现及原理](https://charleswu.site/archives/3065)
- [目前较前卫的 SIMD 优化 NTT](https://www.researchgate.net/publication/358610298_Fast_Implementation_of_Multiplication_on_Polynomial_Rings)
