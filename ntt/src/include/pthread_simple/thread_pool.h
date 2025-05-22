#pragma once
#include <vector>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <future>
#include <atomic>

class ThreadPool
{
public:
    explicit ThreadPool(size_t num_threads);
    ~ThreadPool();

    // 提交一个任务并返回 future
    template <typename Func, typename... Args>
    auto enqueue(Func &&f, Args &&...args) -> std::future<decltype(f(args...))>;

    void wait(); // 等待所有任务完成

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queue_mutex;
    std::condition_variable condition;
    std::atomic<bool> stop;

    std::mutex wait_mutex;
    std::condition_variable wait_condition;
    std::atomic<int> active_tasks;
};

inline ThreadPool::ThreadPool(size_t num_threads) : stop(false), active_tasks(0)
{
    for (size_t i = 0; i < num_threads; ++i)
    {
        workers.emplace_back([this]()
                             {
            while (true) {
                std::function<void()> task;

                {
                    std::unique_lock<std::mutex> lock(queue_mutex);
                    condition.wait(lock, [this]() { return stop || !tasks.empty(); });
                    if (stop && tasks.empty()) return;

                    task = std::move(tasks.front());
                    tasks.pop();
                    ++active_tasks;
                }

                task();

                {
                    std::lock_guard<std::mutex> lock(wait_mutex);
                    if (--active_tasks == 0 && tasks.empty()) {
                        wait_condition.notify_all();
                    }
                }
            } });
    }
}

inline ThreadPool::~ThreadPool()
{
    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for (auto &worker : workers)
    {
        if (worker.joinable())
            worker.join();
    }
}

template <typename Func, typename... Args>
auto ThreadPool::enqueue(Func &&f, Args &&...args) -> std::future<decltype(f(args...))>
{
    using return_type = decltype(f(args...));
    auto task = std::make_shared<std::packaged_task<return_type()>>(
        std::bind(std::forward<Func>(f), std::forward<Args>(args)...));

    {
        std::lock_guard<std::mutex> lock(queue_mutex);
        tasks.emplace([task]()
                      { (*task)(); });
    }

    condition.notify_one();
    return task->get_future();
}

inline void ThreadPool::wait()
{
    std::unique_lock<std::mutex> lock(wait_mutex);
    wait_condition.wait(lock, [this]()
                        { return tasks.empty() && active_tasks == 0; });
}
