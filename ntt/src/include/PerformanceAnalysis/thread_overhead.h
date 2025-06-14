#pragma once

#include <thread>
#include <vector>
#include <chrono>
#include <functional>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include "../config.h"
#include "../Barrett/op.h"

/**
 * @brief 线程池类，用于复用线程避免重复创建开销
 */
class ThreadPool
{
private:
    std::vector<std::thread> workers;
    std::vector<std::function<void()>> tasks;
    std::atomic<int> tasks_assigned{0};
    std::atomic<int> tasks_completed{0};
    std::mutex tasks_mutex;
    std::condition_variable condition;
    bool stop{false};

public:
    ThreadPool(size_t num_threads)
    {
        for (size_t i = 0; i < num_threads; ++i)
        {
            workers.emplace_back([this]
                                 {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->tasks_mutex);
                        this->condition.wait(lock, [this] { 
                            return this->stop || !this->tasks.empty(); 
                        });
                        
                        if (this->stop && this->tasks.empty()) return;
                        
                        if (!this->tasks.empty()) {
                            task = std::move(this->tasks.back());
                            this->tasks.pop_back();
                        }
                    }
                    
                    if (task) {
                        task();
                        tasks_completed.fetch_add(1);
                    }
                } });
        }
    }

    void enqueue_tasks(const std::vector<std::function<void()>> &new_tasks)
    {
        {
            std::unique_lock<std::mutex> lock(tasks_mutex);
            for (const auto &task : new_tasks)
            {
                tasks.emplace_back(task);
            }
            tasks_assigned = new_tasks.size();
            tasks_completed = 0;
        }
        condition.notify_all();

        // 等待所有任务完成
        while (tasks_completed.load() < tasks_assigned.load())
        {
            std::this_thread::yield();
        }
    }

    ~ThreadPool()
    {
        {
            std::unique_lock<std::mutex> lock(tasks_mutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread &worker : workers)
            worker.join();
    }
};

/**
 * @brief 线程创建开销测试结果结构体
 */
struct ThreadOverheadResult
{
    double thread_pool_time_ms;     // 线程池版本耗时(毫秒)
    double thread_creation_time_ms; // 重复创建版本耗时(毫秒)
    double overhead_ratio;          // 开销比例 (creation_time / pool_time)
    int num_threads;                // 线程数
    int num_iterations;             // 迭代次数
};

/**
 * @brief 使用线程池执行NTT计算的函数
 */
template <typename T>
void ntt_with_thread_pool(T *a, T n, T p, T omega, int num_threads)
{
    BarrettMod<T> mod(p, 100);
    ThreadPool pool(num_threads);

    for (T mid = 1; mid < n; mid <<= 1)
    {
        T Wn = mod.pow(omega, (p - 1) / (mid << 1));

        // 准备任务
        std::vector<std::function<void()>> tasks;
        for (T j = 0; j < n; j += (mid << 1))
        {
            tasks.emplace_back([=, &mod]()
                               {
                T w = 1;
                for (T k = 0; k < mid; ++k, w = mod.mul(w, Wn)) {
                    T x = a[j + k];
                    T y = mod.mul(w, a[j + k + mid]);
                    a[j + k] = mod.add(x, y);
                    a[j + k + mid] = mod.sub(x, y);
                } });
        }

        // 执行任务
        pool.enqueue_tasks(tasks);
    }
}

/**
 * @brief 使用重复创建线程执行NTT计算的函数
 */
template <typename T>
void ntt_with_thread_creation(T *a, T n, T p, T omega, int num_threads)
{
    BarrettMod<T> mod(p, 100);

    for (T mid = 1; mid < n; mid <<= 1)
    {
        T Wn = mod.pow(omega, (p - 1) / (mid << 1));

        // 计算任务数
        std::vector<std::pair<T, T>> task_ranges;
        for (T j = 0; j < n; j += (mid << 1))
        {
            task_ranges.emplace_back(j, mid);
        }

        // 分配任务到线程
        int tasks_per_thread = std::max(1, (int)task_ranges.size() / num_threads);
        std::vector<std::thread> threads;

        for (int t = 0; t < num_threads; ++t)
        {
            int start_task = t * tasks_per_thread;
            int end_task = (t == num_threads - 1) ? task_ranges.size() : (t + 1) * tasks_per_thread;

            if (start_task < task_ranges.size())
            {
                threads.emplace_back([=, &mod, &task_ranges]()
                                     {
                    for (int task_idx = start_task; task_idx < end_task; ++task_idx) {
                        T j = task_ranges[task_idx].first;
                        T task_mid = task_ranges[task_idx].second;
                        T w = 1;
                        for (T k = 0; k < task_mid; ++k, w = mod.mul(w, Wn)) {
                            T x = a[j + k];
                            T y = mod.mul(w, a[j + k + mid]);
                            a[j + k] = mod.add(x, y);
                            a[j + k + mid] = mod.sub(x, y);
                        }
                    } });
            }
        }

        // 等待所有线程完成
        for (auto &thread : threads)
        {
            thread.join();
        }
    }
}

/**
 * @brief 测试线程创建开销
 */
template <typename T>
ThreadOverheadResult measure_thread_overhead(T n, T p, T omega, int num_threads, int iterations = 10)
{
    ThreadOverheadResult result;
    result.num_threads = num_threads;
    result.num_iterations = iterations;

    // 准备测试数据
    std::vector<T> test_data_pool(n), test_data_creation(n);
    for (T i = 0; i < n; ++i)
    {
        test_data_pool[i] = test_data_creation[i] = i % p;
    }

    // 测试线程池版本
    {
        TIMER_START();
        for (int i = 0; i < iterations; ++i)
        {
            // 重置数据
            for (T j = 0; j < n; ++j)
            {
                test_data_pool[j] = j % p;
            }
            ntt_with_thread_pool(test_data_pool.data(), n, p, omega, num_threads);
        }
        TIMER_END();
        result.thread_pool_time_ms = TIMER_ELAPSED() / iterations;
    }

    // 测试重复创建版本
    {
        TIMER_START();
        for (int i = 0; i < iterations; ++i)
        {
            // 重置数据
            for (T j = 0; j < n; ++j)
            {
                test_data_creation[j] = j % p;
            }
            ntt_with_thread_creation(test_data_creation.data(), n, p, omega, num_threads);
        }
        TIMER_END();
        result.thread_creation_time_ms = TIMER_ELAPSED() / iterations;
    }

    result.overhead_ratio = result.thread_creation_time_ms / result.thread_pool_time_ms;

    return result;
}