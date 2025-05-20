#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <arm_neon.h>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <pthread.h>
#include <limits.h>
// 可以自行添加需要的头文件
typedef long long LL;

const int ROOT = 3;
const int BITE = 63;
const int NUM_THREADS = 4;
// const int MOD = 998244353;

void fRead(LL *a, LL *b, int *n, LL *p, int input_id)
{
    // 数据输入函数
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strin = str1 + str2 + ".in";
    char data_path[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), data_path);
    data_path[strin.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    fin >> *n >> *p;
    for (int i = 0; i < *n; i++)
    {
        fin >> a[i];
    }
    for (int i = 0; i < *n; i++)
    {
        fin >> b[i];
    }
}

void fCheck(LL *ab, int n, int input_id)
{
    // 判断多项式乘法结果是否正确
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    for (int i = 0; i < n * 2 - 1; i++)
    {
        long long x;
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

void fWrite(LL *ab, int n, int input_id)
{
    // 数据输出函数, 可以用来输出最终结果, 也可用于调试时输出中间数组
    std::string str1 = "files/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++)
    {
        fout << ab[i] << '\n';
    }
}

void poly_multiply(int *a, int *b, int *ab, int n, int p)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            ab[i + j] = (1LL * a[i] * b[j] % p + ab[i + j]) % p;
        }
    }
}

LL qpow(LL a, LL b, LL p)
{
    LL ans = 1;
    while (b)
    {
        if (b & 1)
            ans = (1LL * ans * a) % p;
        a = (1LL * a * a) % p;
        b >>= 1;
    }
    return ans;
}

void bit_reverse(LL *a, int n)
{
    int j = 0;
    for (int i = 1; i < n; i++)
    {
        int bit = n >> 1;
        while (j & bit)
        {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j)
            std::swap(a[i], a[j]);
    }
}

void ntt(LL *a, int n, LL MOD, bool invert = false)
{
    bit_reverse(a, n);

    for (int len = 2; len <= n; len <<= 1)
    {
        LL wn = qpow(ROOT, (MOD - 1) / len, MOD);
        if (invert)
            wn = qpow(wn, MOD - 2, MOD);

        for (int i = 0; i < n; i += len)
        {
            LL w = 1;
            for (int j = 0; j < len / 2; j++)
            {
                LL u = a[i + j], v = 1LL * w * a[i + j + len / 2] % MOD;
                a[i + j] = (u + v) % MOD;
                a[i + j + len / 2] = (u - v + MOD) % MOD;
                w = 1LL * w * wn % MOD;
            }
        }
    }

    if (invert)
    {
        LL inv_n = qpow(n, MOD - 2, MOD);
        for (int i = 0; i < n; i++)
            a[i] = 1LL * a[i] * inv_n % MOD;
    }
}

struct partial_ntt_arg
{
    LL *a;
    int start;
    int end;
    int len;
    LL MOD;
    LL *w;
};

void *partial_ntt(void *arg)
{
    partial_ntt_arg *narg = (partial_ntt_arg *)arg;
    LL *a = narg->a;
    int start = narg->start;
    int end = narg->end;
    int len = narg->len;
    LL MOD = narg->MOD;
    LL *w = narg->w;
    for (int i = start; i < end; i += len)
    {
        for (int j = 0; j < len / 2; j++)
        {
            LL u = a[i + j], v = 1LL * w[j] * a[i + j + len / 2] % MOD;
            a[i + j] = (u + v) % MOD;
            a[i + j + len / 2] = (u - v + MOD) % MOD;
        }
    }
    pthread_exit(NULL);
}

void pthread_ntt(LL *a, int n, LL MOD, bool invert = false)
{
    bit_reverse(a, n);

    for (int len = 2; len <= n; len <<= 1)
    {
        LL wn = qpow(ROOT, (MOD - 1) / len, MOD);
        if (invert)
            wn = qpow(wn, MOD - 2, MOD);
        LL *w = new LL[len / 2];
        w[0] = 1;
        for (int i = 1; i < len / 2; i++)
            w[i] = 1LL * wn * w[i - 1] % MOD;
        partial_ntt_arg arg[NUM_THREADS];
        int partials[NUM_THREADS + 1];
        partials[0] = 0;
        int block;
        int temp_n = n;
        for (int i = 1; i < NUM_THREADS; i++)
        {
            block = temp_n / len;
            partials[i] = partials[i - 1] + block * len;
            temp_n -= block * len;
        }
        partials[NUM_THREADS] = n;
        for (int i = 0; i < NUM_THREADS; i++)
        {
            arg[i].a = a;
            arg[i].start = partials[i];
            arg[i].end = partials[i + 1];
            arg[i].len = len;
            arg[i].MOD = MOD;
            arg[i].w = w;
        }
        void *status;
        pthread_t thread[NUM_THREADS];
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        for (int i = 0; i < NUM_THREADS; i++)
        {
            pthread_create(&thread[i], &attr, partial_ntt, (void *)&arg[i]);
        }
        pthread_attr_destroy(&attr);
        for (int i = 0; i < NUM_THREADS; i++)
        {
            pthread_join(thread[i], &status);
        }
        delete[] w;
    }

    if (invert)
    {
        LL inv_n = qpow(n, MOD - 2, MOD);
        for (int i = 0; i < n; i++)
            a[i] = 1LL * a[i] * inv_n % MOD;
    }
}

void openmp_ntt(LL *a, int n, LL MOD, bool invert = false)
{
    bit_reverse(a, n);

    for (int len = 2; len <= n; len <<= 1)
    {
        LL wn = qpow(ROOT, (MOD - 1) / len, MOD);
        if (invert)
            wn = qpow(wn, MOD - 2, MOD);
        LL *w = new LL[len / 2];
        w[0] = 1;
        for (int i = 1; i < len / 2; i++)
            w[i] = 1LL * wn * w[i - 1] % MOD;

        omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
        for (int i = 0; i < n; i += len)
        {

            for (int j = 0; j < len / 2; j++)
            {
                LL u = a[i + j], v = 1LL * w[j] * a[i + j + len / 2] % MOD;
                a[i + j] = (u + v) % MOD;
                a[i + j + len / 2] = (u - v + MOD) % MOD;
            }
        }
    }

    if (invert)
    {
        LL inv_n = qpow(n, MOD - 2, MOD);
        for (int i = 0; i < n; i++)
            a[i] = 1LL * a[i] * inv_n % MOD;
    }
}
void openmp_ntt_multiply(LL *a, LL *b, LL *ab, int n, LL MOD)
{
    int size = 1;
    while (size < 2 * n)
        size <<= 1;
    for (int i = n; i < size; i++)
        a[i] = b[i] = 0;
    openmp_ntt(a, size, MOD, false);
    openmp_ntt(b, size, MOD, false);
    for (int i = 0; i < size; i++)
        ab[i] = 1LL * a[i] * b[i] % MOD;
    openmp_ntt(ab, size, MOD, true);
}

void pthread_ntt_multiply(LL *a, LL *b, LL *ab, int n, LL MOD)
{
    int size = 1;
    while (size < 2 * n)
        size <<= 1;
    for (int i = n; i < size; i++)
        a[i] = b[i] = 0;
    pthread_ntt(a, size, MOD, false);
    pthread_ntt(b, size, MOD, false);
    for (int i = 0; i < size; i++)
        ab[i] = 1LL * a[i] * b[i] % MOD;
    pthread_ntt(ab, size, MOD, true);
}

void ntt_multiply(LL *a, LL *b, LL *ab, int n, LL MOD)
{
    int size = 1;
    while (size < 2 * n)
        size <<= 1;
    for (int i = n; i < size; i++)
        a[i] = b[i] = 0;
    ntt(a, size, MOD, false);
    ntt(b, size, MOD, false);
    for (int i = 0; i < size; i++)
        ab[i] = 1LL * a[i] * b[i] % MOD;
    ntt(ab, size, MOD, true);
}

struct ntt_arg
{
    LL *a;
    LL *b;
    LL *ab;
    int n;
    LL p; // 和MOD一样
};

void *pthread_ntt_multiply(void *arg)
{
    ntt_arg *narg = (ntt_arg *)arg;
    LL *a = narg->a;
    LL *b = narg->b;
    LL *ab = narg->ab;
    int n = narg->n;
    LL MOD = narg->p;
    int size = 1;
    while (size < 2 * n)
        size <<= 1;
    for (int i = n; i < size; i++)
        a[i] = b[i] = 0;
    ntt(a, size, MOD, false);
    ntt(b, size, MOD, false);
    for (int i = 0; i < size; i++)
        ab[i] = 1LL * a[i] * b[i] % MOD;
    ntt(ab, size, MOD, true);
    pthread_exit(NULL);
}

LL extend_gcd(LL a, LL b, LL &x, LL &y)
{
    if (b == 0)
    {
        x = 1;
        y = 0;
        return a;
    }
    LL x1, y1;
    LL d = extend_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return d;
}

LL inv(LL a, LL p)
{
    LL x, y;
    LL d = extend_gcd(a, p, x, y);
    if (d != 1)
        return -1;
    return (x % p + p) % p;
}

inline LL mulmod(LL a, LL b, LL p)
{
    unsigned long long x = 1ULL * a % p, y = 1ULL * b % p;
    LL result = x * y - (LL)((long double)x / p * y + 0.5) * p;
    return result < 0 ? result + p : result;
    // return 1ULL * a * b % p;
}
// void CRT(LL *ab1, LL *ab2, int n, int p1, int p2){
//     /*
//     * 中国剩余定理合并，最终计算结果会保存在ab1里
//     * @param n为多项式长度
//     * @param p1,p2为两个多项式的模数
//     * @param ab1,ab2为两个多项式的系数
//     */
//     LL P = 1LL * p1 * p2;
//     LL P1 = P / p1;
//     LL P2 = P / p2;
//     LL inv_P1 = inv(P1, (LL)(p1));
//     LL inv_P2 = inv(P2, (LL)(p2));
//     LL temp1 = mulmod(P1 , inv_P1 , P);
//     LL temp2 = mulmod(P2 , inv_P2 , P);
//     LL m1 = temp1 / P1;
//     LL m2 = temp2 / P2;
//     bool flag = false;
//     for(int i = 0; i < n; i++){
//         // if (flag == false && (LONG_LONG_MAX/temp1 < (ab1[i]% P) || LONG_LONG_MAX/temp2 < (ab2[i]% P))){
//         //     std::cout<<"CRT overflow"<<std::endl<<"temp1:"<<temp1<<",temp2:"<<temp2<<",ab1[i]:"<<ab1[i]<<",ab2[i]:"<<ab2[i]<<"p1:"<<p1<<",p2:"<<p2<<std::endl;
//         //     flag = true;
//         // }
//         ab1[i] = ((mulmod(ab1[i], temp1, P) + mulmod(ab2[i], temp2, P))% P+P)%P;
//         // ab1[i] = (1LL * m1 * ab1[i] % P2 * P1 + (1LL * m2 * ab2[i] % P1 * P2)) % P;
//     }
// }
void CRT(LL *ab1, LL *ab2, int n, int p1, int p2)
{
    /*
     * 中国剩余定理合并，最终计算结果会保存在ab1里
     * @param n为多项式长度
     * @param p1,p2为两个多项式的模数
     * @param ab1,ab2为两个多项式的系数
     */
    LL P = 1LL * p1 * p2;
    LL P1 = P / p1;
    LL P2 = P / p2;
    LL inv_P1 = inv(P1, (LL)(p1));
    LL inv_P2 = inv(P2, (LL)(p2));
    LL temp1 = inv_P1 % P2 * P1;
    LL temp2 = inv_P2 % P1 * P2;
    LL m1 = temp1 / P1;
    LL m2 = temp2 / P2;
    bool flag = false;
    for (int i = 0; i < n; i++)
    {
        // if (flag == false && (LONG_LONG_MAX/temp1 < (ab1[i]% P) || LONG_LONG_MAX/temp2 < (ab2[i]% P))){
        //     std::cout<<"CRT overflow"<<std::endl<<"temp1:"<<temp1<<",temp2:"<<temp2<<",ab1[i]:"<<ab1[i]<<",ab2[i]:"<<ab2[i]<<"p1:"<<p1<<",p2:"<<p2<<std::endl;
        //     flag = true;
        // }
        // ab1[i] = ((mulmod(ab1[i], temp1, P) + mulmod(ab2[i], temp2, P))% P+P)%P;
        ab1[i] = (1LL * m1 * ab1[i] % P2 * P1 + (1LL * m2 * ab2[i] % P1 * P2)) % P;
    }
}

void pthread_crt_ntt_multiply(LL *a, LL *b, LL *ab, int n, LL p)
{
    // 计算长度，用于创建复制变量
    int size = 1;
    while (size < 2 * n)
        size <<= 1;
    LL *a_copy[4];
    LL *b_copy[4];
    LL *ab_copy[4];
    struct ntt_arg arg[4];
    // TODO：寻找合适的模数，能找到四个都是3的吗......还真能找到：
    //  https://blog.miskcoo.com/2014/07/fft-prime-table 记录了常用的素数及其原根，致谢@miskcoo
    LL ntt_p[4] = {167772161, 469762049, 998244353, 1004535809};
    a_copy[0] = a;
    b_copy[0] = b;
    ab_copy[0] = ab;
    for (int i = 1; i <= 3; i++)
    {
        a_copy[i] = new LL[size];
        b_copy[i] = new LL[size];
        ab_copy[i] = new LL[size];
    }
    for (int i = 1; i <= 3; i++)
    {
        std::copy(a, a + size, a_copy[i]);
        std::copy(b, b + size, b_copy[i]);
        std::fill(ab_copy[i], ab_copy[i] + size, 0);
    }
    for (int i = 0; i <= 3; i++)
    {
        arg[i].a = a_copy[i];
        arg[i].b = b_copy[i];
        arg[i].ab = ab_copy[i];
        arg[i].n = n;
        arg[i].p = ntt_p[i];
    }
    void *status;
    pthread_t threads[NUM_THREADS];
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    for (int i = 0; i < NUM_THREADS; i++)
    {
        pthread_create(&threads[i], &attr, pthread_ntt_multiply, (void *)&arg[i]);
    }
    pthread_attr_destroy(&attr);
    for (int i = 0; i < NUM_THREADS; i++)
    {
        pthread_join(threads[i], &status);
    }
    for (int i = 0; i < 2; i++)
    {
        CRT(ab_copy[2 * i], ab_copy[2 * i + 1], size, ntt_p[2 * i], ntt_p[2 * i + 1]);
    }

    LL p1 = 1LL * ntt_p[0] * ntt_p[1];
    LL p2 = 1LL * ntt_p[2] * ntt_p[3];
    LL p2_temp = p2 % p;
    LL inv_p2 = inv(p2, p1);
    for (int i = 0; i < size; i++)
    {
        // LL k = ((ab_copy[0][i]-ab_copy[2][i])% p1 + p1)%p1 * (inv_p2 % p1) % p1;
        LL k = mulmod(((ab_copy[0][i] - ab_copy[2][i]) % p1 + p1) % p1, (inv_p2 % p1), p1);
        LL temp = (mulmod(k, p2, p) + ab_copy[2][i]) % p;
        ab[i] = temp;
    }
    for (int i = 1; i <= 3; i++)
    {
        delete[] a_copy[i];
        delete[] b_copy[i];
        delete[] ab_copy[i];
    }
    // for(int i=0;i<1;i++){
    //     CRT(ab_copy[2*i], ab_copy[2*i+1], size, ntt_p[2*i], ntt_p[2*i+1]);
    // }
    // LL p1 = 1LL * ntt_p[0] * ntt_p[1];
    // LL p2 = 1LL * ntt_p[2];
    // LL p1_temp = p1 % p;
    // LL inv_p1 = inv(p1, p2);
    // for(int i=0;i<size;i++){
    //     LL k = 1LL * (((ab_copy[2][i]-ab_copy[0][i])% p2 + p2)%p2)* (inv_p1 % p2) % p2;
    //     ab[i] = ((((k%p)*p1_temp)%p + (LL)ab_copy[0][i]%p)%p+p)%p;
    // }
    // for(int i=1;i<=3;i++){
    //     delete[] a_copy[i];
    //     delete[] b_copy[i];
    //     delete[] ab_copy[i];
    // }
}

void openmp_crt_ntt_multiply(LL *a, LL *b, LL *ab, int n, LL p)
{
    // 计算长度，用于创建复制变量
    int size = 1;
    while (size < 2 * n)
        size <<= 1;
    LL *a_copy[4];
    LL *b_copy[4];
    LL *ab_copy[4];
    // TODO：寻找合适的模数，能找到四个都是3的吗......还真能找到：
    //  https://blog.miskcoo.com/2014/07/fft-prime-table 记录了常用的素数及其原根，致谢@miskcoo
    LL ntt_p[4] = {167772161, 469762049, 998244353, 1004535809};
    a_copy[0] = a;
    b_copy[0] = b;
    ab_copy[0] = ab;
    for (int i = 1; i <= 3; i++)
    {
        a_copy[i] = new LL[size];
        b_copy[i] = new LL[size];
        ab_copy[i] = new LL[size];
    }
    for (int i = 1; i <= 3; i++)
    {
        std::copy(a, a + size, a_copy[i]);
        std::copy(b, b + size, b_copy[i]);
        std::fill(ab_copy[i], ab_copy[i] + size, 0);
    }
    omp_set_num_threads(NUM_THREADS);
#pragma omp parallel for
    for (int i = 0; i <= 3; i++)
    {
        ntt_multiply(a_copy[i], b_copy[i], ab_copy[i], n, ntt_p[i]);
    }
    for (int i = 0; i < 2; i++)
    {
        CRT(ab_copy[2 * i], ab_copy[2 * i + 1], size, ntt_p[2 * i], ntt_p[2 * i + 1]);
    }

    LL p1 = 1LL * ntt_p[0] * ntt_p[1];
    LL p2 = 1LL * ntt_p[2] * ntt_p[3];
    LL p2_temp = p2 % p;
    LL inv_p2 = inv(p2, p1);
    // omp_set_num_threads(NUM_THREADS);
    // #pragma omp parallel for
    for (int i = 0; i < size; i++)
    {
        // LL k = ((ab_copy[0][i]-ab_copy[2][i])% p1 + p1)%p1 * (inv_p2 % p1) % p1;
        LL k = mulmod(((ab_copy[0][i] - ab_copy[2][i]) % p1 + p1) % p1, (inv_p2 % p1), p1);
        LL temp = (mulmod(k, p2, p) + ab_copy[2][i]) % p;
        ab[i] = temp;
    }
    for (int i = 1; i <= 3; i++)
    {
        delete[] a_copy[i];
        delete[] b_copy[i];
        delete[] ab_copy[i];
    }
}

void crt_ntt_multiply(LL *a, LL *b, LL *ab, int n, LL p)
{
    // 计算长度，用于创建复制变量
    int size = 1;
    while (size < 2 * n)
        size <<= 1;
    LL *a_copy[4];
    LL *b_copy[4];
    LL *ab_copy[4];
    // TODO：寻找合适的模数，能找到四个都是3的吗......还真能找到：
    //  https://blog.miskcoo.com/2014/07/fft-prime-table 记录了常用的素数及其原根，致谢@miskcoo
    LL ntt_p[4] = {167772161, 469762049, 998244353, 1004535809};
    a_copy[0] = a;
    b_copy[0] = b;
    ab_copy[0] = ab;
    for (int i = 1; i <= 3; i++)
    {
        a_copy[i] = new LL[size];
        b_copy[i] = new LL[size];
        ab_copy[i] = new LL[size];
    }
    for (int i = 1; i <= 3; i++)
    {
        std::copy(a, a + size, a_copy[i]);
        std::copy(b, b + size, b_copy[i]);
        std::fill(ab_copy[i], ab_copy[i] + size, 0);
    }
    for (int i = 0; i <= 3; i++)
    {
        ntt_multiply(a_copy[i], b_copy[i], ab_copy[i], n, ntt_p[i]);
    }
    for (int i = 0; i < 2; i++)
    {
        CRT(ab_copy[2 * i], ab_copy[2 * i + 1], size, ntt_p[2 * i], ntt_p[2 * i + 1]);
    }

    LL p1 = 1LL * ntt_p[0] * ntt_p[1];
    LL p2 = 1LL * ntt_p[2] * ntt_p[3];
    LL p2_temp = p2 % p;
    LL inv_p2 = inv(p2, p1);
    for (int i = 0; i < size; i++)
    {
        // LL k = ((ab_copy[0][i]-ab_copy[2][i])% p1 + p1)%p1 * (inv_p2 % p1) % p1;
        LL k = mulmod(((ab_copy[0][i] - ab_copy[2][i]) % p1 + p1) % p1, (inv_p2 % p1), p1);
        LL temp = (mulmod(k, p2, p) + ab_copy[2][i]) % p;
        ab[i] = temp;
    }
    for (int i = 1; i <= 3; i++)
    {
        delete[] a_copy[i];
        delete[] b_copy[i];
        delete[] ab_copy[i];
    }
}

inline void neon_mod(int32x4_t *a, int64_t mu, int MOD)
{
    int32_t vals[4];
    int32_t subs[4];
    vst1q_s32(vals, *a);
    for (int i = 0; i < 4; i++)
    {
        int64_t q = (vals[i] * mu) >> BITE;
        // int64_t q = vals[i]/ MOD;
        subs[i] = static_cast<int32_t>(q * MOD);
    }
    int32x4_t sub = vld1q_s32(subs);
    *a = vsubq_s32(*a, sub);

    int i = 0;
    while (true)
    {
        int32x4_t mod = vdupq_n_s32((1 << i) * MOD);
        uint32x4_t mask = vcgeq_s32(*a, mod);
        uint32_t flag = vmaxvq_u32(mask);
        if (flag == 0)
            break;
        ++i;
    }
    if (i > 0)
        i--;
    while (true)
    {
        int32x4_t mod = vdupq_n_s32((1 << i) * MOD);
        uint32x4_t mask = vcgeq_s32(*a, mod);
        uint32_t flag = vmaxvq_u32(mask);
        if (flag)
        {
            int32x4_t subtracted = vsubq_s32(*a, mod);
            *a = vbslq_s32(mask, subtracted, *a);
        }
        if (i == 0)
            break;
        --i;
    }
}

inline uint64_t vmaxvq_u64(uint64x2_t v)
{
    uint64_t a = vgetq_lane_u64(v, 0);
    uint64_t b = vgetq_lane_u64(v, 1);
    return a > b ? a : b;
}

uint64x2_t neon_mul_int64(int64_t a, int64_t b)
{
    a = uint64_t(a);
    b = uint64_t(b);

    uint32_t a_lo = (uint32_t)a;
    uint32_t a_hi = (uint32_t)(a >> 32);
    uint32_t b_lo = (uint32_t)b;
    uint32_t b_hi = (uint32_t)(b >> 32);

    uint64x2_t result = vdupq_n_u64(0);

    uint64_t low = (uint64_t)a_lo * b_lo;

    uint64_t mid = (uint64_t)a_lo * b_hi + (uint64_t)a_hi * b_lo;
    uint32_t mid_hi = (uint32_t)(mid >> 32);
    uint32_t mid_lo = (uint32_t)mid;

    uint64_t high = (uint64_t)a_hi * b_hi;

    uint64_t lo = (uint64_t)mid_lo << 32;
    uint64_t final_low = low + lo;

    uint64_t carry = (0xFFFFFFFFFFFFFFFF - lo < low);

    uint64_t final_high = high + (mid_hi) + carry;

    result = vsetq_lane_u64(final_low, result, 0);
    result = vsetq_lane_u64(final_high, result, 1);
    return result;
}

inline void neon_mod_for_mul(int64x2_t *a, int MOD)
{
    int64x2_t mod = vdupq_n_s64(MOD);
    uint64x2_t mask = vcgeq_s64(*a, mod);
    int64x2_t subtracted = vsubq_s64(*a, mod);
    *a = vbslq_s64(mask, subtracted, *a);
}

inline void neon_barrett_reduction(int64x2_t *a, int64_t mu, int MOD)
{
    int64_t vals[2];
    int64_t subs[2];
    vst1q_s64(vals, *a);
    for (int i = 0; i < 2; i++)
    {
        uint64x2_t q_temp = neon_mul_int64(vals[i], mu);
        uint64_t q_lo = vgetq_lane_u64(q_temp, 0);
        uint64_t q_hi = vgetq_lane_u64(q_temp, 1);
        int64_t q = (q_lo >> BITE) + (q_hi << (64 - BITE));
        subs[i] = q * MOD;
    }
    int64x2_t sub = vld1q_s64(subs);
    *a = vsubq_s64(*a, sub);
    int64x2_t mod = vdupq_n_s64(MOD);
    uint64x2_t mask = vcgeq_s64(*a, mod);
    int64x2_t subtracted = vsubq_s64(*a, mod);
    *a = vbslq_s64(mask, subtracted, *a);
}

inline int32x4_t neon_add(int32x4_t *a, int32x4_t *b, int MOD, int64_t MU)
{
    int32x4_t c = vaddq_s32(*a, *b);
    neon_mod(&c, MU, MOD);
    return c;
}

inline int32x4_t neon_sub(int32x4_t *a, int32x4_t *b, int MOD, int64_t MU)
{
    int32x4_t c = vsubq_s32(*a, *b);
    c = vaddq_s32(c, vdupq_n_s32(MOD));
    neon_mod(&c, MU, MOD);
    return c;
}

inline int32x4_t neon_mul(int32x4_t *a, int32x4_t *b, int MOD, int64_t MU)
{
    // int32x4_t c = vmulq_s32(*a, *b);
    // neon_mod(&c, MOD);
    // return c;
    int64x2_t c_lo = vmull_s32(vget_low_s32(*a), vget_low_s32(*b));
    int64x2_t c_hi = vmull_s32(vget_high_s32(*a), vget_high_s32(*b));
    neon_barrett_reduction(&c_lo, MU, MOD);
    neon_barrett_reduction(&c_hi, MU, MOD);
    return vcombine_s32(vmovn_s64(c_lo), vmovn_s64(c_hi));
}

void neon_ntt(int *a, int n, int MOD, int64_t MU, bool invert = false)
{
    // bit_reverse(a, n);

    for (int len = 2; len <= n; len <<= 1)
    {
        int wn = qpow(ROOT, (MOD - 1) / len, MOD);
        if (invert)
            wn = qpow(wn, MOD - 2, MOD);

        for (int i = 0; i < n; i += len)
        {
            int w = 1;
            for (int j = 0; j < len / 2; j += 4)
            {
                if (j + 3 < len / 2)
                {
                    int32_t temp[4];
                    for (int k = 0; k < 4; k++)
                    {
                        temp[k] = w;
                        w = 1LL * w * wn % MOD;
                    }
                    int32x4_t ww = vld1q_s32(temp);
                    int32x4_t u = vld1q_s32(a + i + j);
                    int32x4_t v = vld1q_s32(a + i + j + len / 2);
                    int32x4_t wv = neon_mul(&ww, &v, MOD, MU);

                    int32x4_t u_add_v = neon_add(&u, &wv, MOD, MU);
                    int32x4_t u_sub_v = neon_sub(&u, &wv, MOD, MU);

                    vst1q_s32(a + i + j, u_add_v);
                    vst1q_s32(a + i + j + len / 2, u_sub_v);
                }
                else
                {
                    for (int k = j; k < len / 2; k++)
                    {
                        int u = a[i + k], v = 1LL * w * a[i + k + len / 2] % MOD;
                        a[i + k] = (u + v) % MOD;
                        a[i + k + len / 2] = (u - v + MOD) % MOD;
                        w = 1LL * w * wn % MOD;
                    }
                }
            }
        }
    }

    if (invert)
    {
        int inv_n = qpow(n, MOD - 2, MOD);
        for (int i = 0; i < n; i += 4)
        {
            int32x4_t a4 = vld1q_s32(a + i);
            int32x4_t inv_n4 = vdupq_n_s32(inv_n);
            int32x4_t ab4 = neon_mul(&a4, &inv_n4, MOD, MU);
            vst1q_s32(a + i, ab4);
            //   a[i] = 1LL * a[i] * inv_n % MOD;
        }
    }
}

void neon_ntt_multiply(int *a, int *b, int *ab, int n, int MOD)
{
    const int64_t MU = (1LL << BITE) / static_cast<int64_t>(MOD);
    int size = 1;
    while (size < 2 * n)
        size <<= 1;
    for (int i = n; i < size; i++)
        a[i] = b[i] = 0;
    neon_ntt(a, size, MOD, MU, false);
    neon_ntt(b, size, MOD, MU, false);
    for (int i = 0; i < size; i += 4)
    {
        // ab[i] = 1LL * a[i] * b[i] % MOD;
        int32x4_t a4 = vld1q_s32(a + i);
        int32x4_t b4 = vld1q_s32(b + i);
        int32x4_t ab4 = neon_mul(&a4, &b4, MOD, MU);
        vst1q_s32(ab + i, ab4);
    }
    neon_ntt(ab, size, MOD, MU, true);
}

LL a[300000], b[300000], ab[300000];
int main(int argc, char *argv[])
{

    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 1337006139375617
    // 第四个模数超过了整型表示范围, 如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求, 如果要自行探索大模数 NTT, 请在完成前三个模数的基础代码及优化后实现大模数 NTT
    // 输入文件共五个, 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久, 推荐调试正确性时只使用输入文件 1
    int test_begin = 0;
    int test_end = 4;
    for (int i = test_begin; i <= test_end; ++i)
    {
        long double ans = 0;
        int n_;
        LL p_;
        fRead(a, b, &n_, &p_, i);
        memset(ab, 0, sizeof(ab));
        auto Start = std::chrono::high_resolution_clock::now();
        // TODO : 将 poly_multiply 函数替换成你写的 ntt
        // poly_multiply(a, b, ab, n_, p_);
        // ntt_multiply(a, b, ab, n_, p_);
        // neon_ntt_multiply(a, b, ab, n_, p_);
        pthread_crt_ntt_multiply(a, b, ab, n_, p_);
        // openmp_crt_ntt_multiply(a, b, ab, n_, p_);
        // openmp_ntt_multiply(a, b, ab, n_, p_);
        // pthread_ntt_multiply(a, b, ab, n_, p_);
        // crt_ntt_multiply(a, b, ab, n_, p_);
        auto End = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::ratio<1, 1000>> elapsed = End - Start;
        ans += elapsed.count();
        fCheck(ab, n_, i);
        std::cout << "average latency for n = " << n_ << " p = " << p_ << " : " << ans << " (us) " << std::endl;
        // 可以使用 fWrite 函数将 ab 的输出结果打印到 files 文件夹下
        // 禁止使用 cout 一次性输出大量文件内容
        fWrite(ab, n_, i);
    }

    pthread_exit(NULL);
}
