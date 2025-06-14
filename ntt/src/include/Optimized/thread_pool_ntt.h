#pragma once

#include "../config.h"
#include "../Barrett/op.h"
#include <thread>
#include <vector>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <future>
#include <atomic>
#include <iostream>

/**
 * 高性能线程池实现
 * 专门为NTT算法优化，避免频繁的线程创建和销毁
 */
class NTTThreadPool
{
private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queue_mutex;
    std::condition_variable condition;
    std::atomic<bool> stop{false};
    std::atomic<size_t> active_tasks{0};

public:
    explicit NTTThreadPool(size_t threads = std::thread::hardware_concurrency())
    {
        for (size_t i = 0; i < threads; ++i)
        {
            workers.emplace_back([this]
                                 {
                while(true) {
                    std::function<void()> task;
                    
                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock, [this] { 
                            return this->stop || !this->tasks.empty(); 
                        });
                        
                        if(this->stop && this->tasks.empty())
                            return;
                            
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    
                    active_tasks++;
                    task();
                    active_tasks--;
                } });
        }
    }

    template <class F, class... Args>
    auto enqueue(F &&f, Args &&...args) -> std::future<typename std::result_of<F(Args...)>::type>
    {
        using return_type = typename std::result_of<F(Args...)>::type;

        auto task = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<F>(f), std::forward<Args>(args)...));

        std::future<return_type> res = task->get_future();

        {
            std::unique_lock<std::mutex> lock(queue_mutex);
            if (stop)
                throw std::runtime_error("enqueue on stopped ThreadPool");

            tasks.emplace([task]()
                          { (*task)(); });
        }

        condition.notify_one();
        return res;
    }

    void wait_all()
    {
        while (active_tasks > 0 || !tasks.empty())
        {
            std::this_thread::sleep_for(std::chrono::microseconds(10));
        }
    }

    size_t get_thread_count() const
    {
        return workers.size();
    }

    ~NTTThreadPool()
    {
        stop = true;
        condition.notify_all();
        for (std::thread &worker : workers)
            worker.join();
    }
};

// 全局线程池实例
static NTTThreadPool g_ntt_thread_pool;

/**
 * 基于线程池的NTT正变换
 * 消除线程创建开销，复用已创建的线程
 */
template <typename T>
inline void ntt_forward_thread_pool(T *a, T n, T p, T omega)
{
    BarrettMod<T> mod(p, 100);

    for (T mid = 1; mid < n; mid <<= 1)
    {
        T Wn = mod.pow(omega, (p - 1) / (mid << 1));

        // 计算任务数量和每个任务的工作量
        T total_blocks = n / (mid << 1);
        T tasks_per_thread = std::max(T(1), total_blocks / g_ntt_thread_pool.get_thread_count());

        std::vector<std::future<void>> futures;

        // 将工作分配给线程池
        for (T start = 0; start < total_blocks; start += tasks_per_thread)
        {
            T end = std::min(start + tasks_per_thread, total_blocks);

            auto future = g_ntt_thread_pool.enqueue([=, &mod]()
                                                    {
                for (T block = start; block < end; ++block) {
                    T j = block * (mid << 1);
                    T w = 1;
                    
                    for (T k = 0; k < mid; ++k, w = mod.mul(w, Wn))
                    {
                        T x = a[j + k];
                        T y = mod.mul(w, a[j + k + mid]);
                        a[j + k] = mod.add(x, y);
                        a[j + k + mid] = mod.sub(x, y);
                    }
                } });

            futures.push_back(std::move(future));
        }

        // 等待所有任务完成
        for (auto &future : futures)
        {
            future.wait();
        }
    }
}

/**
 * 基于线程池的NTT逆变换
 */
template <typename T>
inline void ntt_backward_thread_pool(T *a, T n, T p, T omega)
{
    BarrettMod<T> mod(p, 100);
    T omega_inv = mod.pow(omega, p - 2);

    for (T mid = n >> 1; mid >= 1; mid >>= 1)
    {
        T Wn = mod.pow(omega_inv, (p - 1) / (mid << 1));

        T total_blocks = n / (mid << 1);
        T tasks_per_thread = std::max(T(1), total_blocks / g_ntt_thread_pool.get_thread_count());

        std::vector<std::future<void>> futures;

        for (T start = 0; start < total_blocks; start += tasks_per_thread)
        {
            T end = std::min(start + tasks_per_thread, total_blocks);

            auto future = g_ntt_thread_pool.enqueue([=, &mod]()
                                                    {
                for (T block = start; block < end; ++block) {
                    T j = block * (mid << 1);
                    T w = 1;
                    
                    for (T k = 0; k < mid; ++k, w = mod.mul(w, Wn))
                    {
                        T x = a[j + k];
                        T y = a[j + k + mid];
                        a[j + k] = mod.add(x, y);
                        a[j + k + mid] = mod.mul(mod.sub(x, y), w);
                    }
                } });

            futures.push_back(std::move(future));
        }

        for (auto &future : futures)
        {
            future.wait();
        }
    }

    // 归一化处理
    T n_inv = mod.pow(n, p - 2);
    T total_elements = n;
    T elements_per_task = std::max(T(1), total_elements / g_ntt_thread_pool.get_thread_count());

    std::vector<std::future<void>> normalize_futures;

    for (T start = 0; start < total_elements; start += elements_per_task)
    {
        T end = std::min(start + elements_per_task, total_elements);

        auto future = g_ntt_thread_pool.enqueue([=, &mod]()
                                                {
            for (T i = start; i < end; ++i) {
                a[i] = mod.mul(a[i], n_inv);
            } });

        normalize_futures.push_back(std::move(future));
    }

    for (auto &future : normalize_futures)
    {
        future.wait();
    }
}

/**
 * 基于线程池的多项式乘法
 */
template <typename T>
inline void poly_multiply_ntt_thread_pool(T *a, T *b, T *ab, T n, T p, T omega)
{
    T n_expanded = expand_n(2 * n - 1);

    T *a_expanded = new T[n_expanded]{};
    T *b_expanded = new T[n_expanded]{};

    // 并行复制数据
    auto copy_future_a = g_ntt_thread_pool.enqueue([=]()
                                                   { std::copy(a, a + n, a_expanded); });

    auto copy_future_b = g_ntt_thread_pool.enqueue([=]()
                                                   { std::copy(b, b + n, b_expanded); });

    copy_future_a.wait();
    copy_future_b.wait();

    // 并行执行NTT正变换
    auto ntt_future_a = g_ntt_thread_pool.enqueue([=]()
                                                  { ntt_forward_thread_pool(a_expanded, n_expanded, p, omega); });

    auto ntt_future_b = g_ntt_thread_pool.enqueue([=]()
                                                  { ntt_forward_thread_pool(b_expanded, n_expanded, p, omega); });

    ntt_future_a.wait();
    ntt_future_b.wait();

    // 点乘
    BarrettMod<T> mod(p);
    T elements_per_task = std::max(T(1), n_expanded / g_ntt_thread_pool.get_thread_count());

    std::vector<std::future<void>> mul_futures;

    for (T start = 0; start < n_expanded; start += elements_per_task)
    {
        T end = std::min(start + elements_per_task, n_expanded);

        auto future = g_ntt_thread_pool.enqueue([=, &mod]()
                                                {
            for (T i = start; i < end; ++i) {
                a_expanded[i] = mod.mul(a_expanded[i], b_expanded[i]);
            } });

        mul_futures.push_back(std::move(future));
    }

    for (auto &future : mul_futures)
    {
        future.wait();
    }

    // 逆变换
    ntt_backward_thread_pool(a_expanded, n_expanded, p, omega);

    // 复制结果
    std::copy(a_expanded, a_expanded + n_expanded, ab);

    delete[] a_expanded;
    delete[] b_expanded;
}

/**
 * 线程池性能统计
 */
struct ThreadPoolStats
{
    std::atomic<size_t> tasks_executed{0};
    std::atomic<size_t> total_wait_time_us{0};
    std::chrono::high_resolution_clock::time_point start_time;

    void reset()
    {
        tasks_executed = 0;
        total_wait_time_us = 0;
        start_time = std::chrono::high_resolution_clock::now();
    }

    void record_task()
    {
        tasks_executed++;
    }

    void record_wait_time(std::chrono::microseconds wait_time)
    {
        total_wait_time_us += wait_time.count();
    }

    void print() const
    {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

        std::cout << "[线程池性能统计]" << std::endl;
        std::cout << "执行任务数: " << tasks_executed.load() << std::endl;
        std::cout << "总等待时间: " << total_wait_time_us.load() << " 微秒" << std::endl;
        std::cout << "总运行时间: " << total_time.count() << " 微秒" << std::endl;
        std::cout << "线程利用率: " << (1.0 - double(total_wait_time_us.load()) / total_time.count()) * 100 << "%" << std::endl;
    }
};

// 全局统计对象
extern ThreadPoolStats g_thread_pool_stats;