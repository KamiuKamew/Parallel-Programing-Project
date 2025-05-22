#pragma once

#include "../general/op.h"
#include "thread_pool.h"
#include <thread>
#include <algorithm>

/**
 * @brief NTT 正变换：a(x) → A(ω)
 *
 * 输入的顺序是 bit-reversed，输出的顺序是自然顺序。
 *
 * 进行就地变换是缓存友好的。
 *
 * @param a_mont 多项式系数（位于 Montgomery 数域），变换后表示频域系数（位于 Montgomery 数域）
 * @param n 多项式长度（普通整数）
 * @param p 模数（普通整数）
 * @param omega_mont 原根（位于 Montgomery 数域）
 */
template <typename T>
inline void ntt_forward_mont_pthread_simple(T *a_mont, T n, T p, T omega_mont)
{
    using T_mont = T;
    MontMod<T> montMod(p);
    size_t num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0)
        num_threads = 1;

    ThreadPool pool(num_threads);

    for (T mid = 1; mid < n; mid <<= 1)
    {
        T_mont Wn_mont = montMod.pow(omega_mont, (p - 1) / (mid << 1));
        // for (T j = 0; j < n; j += (mid << 1))
        // {
        //     auto k_loop = [=, &montMod]()
        //     {
        //         T_mont w_loc_mont = montMod.from_T(1);
        //         for (T k = 0; k < mid; ++k, w_loc_mont = montMod.mul(w_loc_mont, Wn_mont))
        //         {
        //             T_mont x_mont_val = a_mont[j + k];
        //             T_mont y_mont_val = montMod.mul(w_loc_mont, a_mont[j + k + mid]);
        //             a_mont[j + k] = montMod.add(x_mont_val, y_mont_val);
        //             a_mont[j + k + mid] = montMod.sub(x_mont_val, y_mont_val);
        //         }
        //     };
        //     if (mid == 1)
        //         k_loop();
        //     else
        //         pool.enqueue(k_loop);
        // }
        T total_blocks = n / (mid << 1);

        for (size_t t = 0; t < num_threads; ++t)
        {
            T start_block = total_blocks * t / num_threads;
            T end_block = total_blocks * (t + 1) / num_threads;

            pool.enqueue([=, &montMod]()
                         {
                for (T b = start_block; b < end_block; ++b)
                {
                    T j = b * (mid << 1);
                    T_mont w_mont = montMod.from_T(1);
                    for (T k = 0; k < mid; ++k, w_mont = montMod.mul(w_mont, Wn_mont))
                    {
                        T_mont x = a_mont[j + k];
                        T_mont y = montMod.mul(w_mont, a_mont[j + k + mid]);
                        a_mont[j + k] = montMod.add(x, y);
                        a_mont[j + k + mid] = montMod.sub(x, y);
                    }
                } });
        }

        pool.wait();
    }
}

/**
 * @brief NTT 逆变换：A(ω) → a_mont(x)
 *
 * 输入的顺序是自然顺序，输出的顺序是 bit-reversed。
 *
 * @param a_mont 频域系数，变换后表示多项式系数
 * @param n 多项式长度
 * @param p 模数
 * @param inv_omega_val_mont 原根，已经是正变换的 ω，在调用时传 montMod.inv(omega_mont)
 */
template <typename T>
inline void ntt_inverse_mont_pthread_simple(T *a_mont, T n, T p, T inv_omega_val_mont)
{
    using T_mont = T;
    MontMod<T> montMod(p);
    size_t num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0)
        num_threads = 1;

    ThreadPool pool(num_threads);

    for (T mid = n >> 1; mid > 0; mid >>= 1)
    {
        T_mont Wn_mont = montMod.pow(inv_omega_val_mont, (p - 1) / (mid << 1));
        // for (T j = 0; j < n; j += (mid << 1))
        // {
        //     auto k_loop = [=, &montMod]()
        //     {
        //         T_mont w_loc_mont = montMod.from_T(1);
        //         for (T k = 0; k < mid; ++k, w_loc_mont = montMod.mul(w_loc_mont, Wn_mont))
        //         {
        //             T_mont x_mont_val = a_mont[j + k];
        //             T_mont y_mont_val = a_mont[j + k + mid];
        //             a_mont[j + k] = montMod.add(x_mont_val, y_mont_val);
        //             a_mont[j + k + mid] = montMod.mul(w_loc_mont, montMod.sub(x_mont_val, y_mont_val));
        //         }
        //     };
        //     if (mid == 1)
        //         k_loop();
        //     else
        //         pool.enqueue(k_loop);
        // }
        T total_blocks = n / (mid << 1);

        for (size_t t = 0; t < num_threads; ++t)
        {
            T start_block = total_blocks * t / num_threads;
            T end_block = total_blocks * (t + 1) / num_threads;

            pool.enqueue([=, &montMod]()
                         {
                for (T b = start_block; b < end_block; ++b)
                {
                    T j = b * (mid << 1);
                    T_mont w_mont = montMod.from_T(1);
                    for (T k = 0; k < mid; ++k, w_mont = montMod.mul(w_mont, Wn_mont))
                    {
                        T_mont x = a_mont[j + k];
                        T_mont y = a_mont[j + k + mid];
                        a_mont[j + k] = montMod.add(x, y);
                        a_mont[j + k + mid] = montMod.mul(w_mont, montMod.sub(x, y));
                    }
                } });
        }

        pool.wait();
    }

    T_mont inv_n = montMod.inv(montMod.from_T(n));
    for (T i = 0; i < n; ++i)
        a_mont[i] = montMod.mul(a_mont[i], inv_n);
}
