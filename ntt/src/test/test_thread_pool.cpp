#include "../include/pthread_simple/thread_pool.h" // Assuming include paths are set up in the build system
#include <iostream>
#include <vector>
#include <atomic>
#include <thread>
#include <chrono>
#include <mutex>
#include <numeric>
#include <algorithm>

// Global mutex to protect the futures vector
std::mutex futures_mutex;

// Simple function to increment a counter
void increment_counter(std::atomic<int> &counter)
{
    // Simulate some work
    // std::this_thread::sleep_for(std::chrono::milliseconds(1)); // Temporarily commented out for linter/build speed
    counter.fetch_add(1, std::memory_order_relaxed);
}

// Function that throws an exception
void task_that_throws()
{
    // Simulate some work
    // std::this_thread::sleep_for(std::chrono::milliseconds(1)); // Temporarily commented out for linter/build speed
    throw std::runtime_error("Simulated task exception");
}

int main()
{
    const int num_threads_in_pool = 4;
    const int num_tasks = 1000;               // Total number of tasks
    const int num_concurrent_submitters = 10; // Number of threads concurrently submitting tasks

    ThreadPool pool(num_threads_in_pool);
    std::atomic<int> counter(0);
    std::vector<std::future<void>> futures;
    futures.reserve(num_tasks);                       // Pre-allocate memory for the futures vector
    std::atomic<int> submitted_tasks_atomic_count(0); // Atomic counter for actually submitted tasks
    std::atomic<int> completed_tasks_count(0);        // Atomic counter for processed futures
    std::atomic<int> exception_from_task_count(0);    // Atomic counter for exceptions caught from tasks

    std::cout << "Starting ThreadPool safety test..." << std::endl;
    std::cout << "Thread pool size: " << num_threads_in_pool << std::endl;
    std::cout << "Total tasks to submit: " << num_tasks << std::endl;
    std::cout << "Number of concurrent submitter threads: " << num_concurrent_submitters << std::endl;

    std::vector<std::thread> submitters;
    std::atomic<int> task_submission_index(0); // Used to determine task type (normal/throwing)

    for (int i = 0; i < num_concurrent_submitters; ++i)
    {
        submitters.emplace_back([&, submitter_id = i]()
                                {
            while(true) {
                int current_task_idx = task_submission_index.fetch_add(1, std::memory_order_relaxed);
                if (current_task_idx >= num_tasks) {
                    break; // All tasks have been claimed for submission
                }

                try {
                    std::future<void> f;
                    if (current_task_idx % 10 == 0) { // Every 10th task (by submission order) throws an exception
                        f = pool.enqueue(task_that_throws);
                    } else {
                        f = pool.enqueue(increment_counter, std::ref(counter));
                    }
                    
                    {
                        std::lock_guard<std::mutex> lock(futures_mutex); // Protect push_back
                        futures.push_back(std::move(f));
                    }
                    submitted_tasks_atomic_count.fetch_add(1, std::memory_order_relaxed);

                } catch (const std::exception& e) {
                    // This catch block is for exceptions from pool.enqueue itself (e.g., if pool is stopped)
                    std::cerr << "Exception during task submission by submitter " << submitter_id 
                              << " for task index " << current_task_idx << ": " << e.what() << std::endl;
                }
            } });
    }

    // Wait for all submitter threads to finish
    for (std::thread &t : submitters)
    {
        if (t.joinable())
        {
            t.join();
        }
    }

    std::cout << "All " << submitted_tasks_atomic_count.load() << " tasks have been enqueued. Processing futures..." << std::endl;
    if (futures.size() != submitted_tasks_atomic_count.load())
    {
        // This could happen if enqueue itself threw and was caught, preventing a future from being added
        std::cout << "Warning: Futures vector size (" << futures.size()
                  << ") does not match submitted tasks count (" << submitted_tasks_atomic_count.load() << ")." << std::endl;
    }

    // Wait for all tasks to complete by getting their futures, and handle exceptions
    for (size_t k = 0; k < futures.size(); ++k)
    {
        try
        {
            if (futures[k].valid())
            {
                futures[k].get(); // Wait for the task to complete
            }
            else
            {
                // This case should ideally not happen if futures are added correctly and only processed once.
                // std::cout << "Warning: Future at index " << k << " is not valid before get()." << std::endl;
            }
        }
        catch (const std::runtime_error &e)
        {
            // This is the expected exception from task_that_throws
            exception_from_task_count.fetch_add(1, std::memory_order_relaxed);
        }
        catch (const std::future_error &e)
        {
            std::cerr << "Future error for future at index " << k << ": " << e.what() << " (code: " << e.code() << ")" << std::endl;
            // This might happen if future.get() is called more than once, or the future is otherwise invalid.
        }
        catch (const std::exception &e)
        {
            std::cerr << "Future at index " << k << " threw an unexpected non-runtime_error exception: " << e.what() << std::endl;
        }
        completed_tasks_count.fetch_add(1, std::memory_order_relaxed); // Increment for each future processed
    }

    std::cout << "All futures processed. Calling pool.wait() to ensure pool is empty..." << std::endl;
    pool.wait(); // Wait for any remaining tasks in the pool (should be none if all futures were processed)
    std::cout << "pool.wait() finished." << std::endl;

    int actual_submitted_tasks = submitted_tasks_atomic_count.load();
    int actual_processed_futures = completed_tasks_count.load();
    int actual_counter_value = counter.load();
    int actual_exceptions_from_tasks = exception_from_task_count.load();

    int expected_normal_tasks = 0;
    int expected_throwing_tasks = 0;

    // Calculate expected values based on the submission logic
    // This assumes tasks were indexed 0 to actual_submitted_tasks-1 when deciding type
    for (int k = 0; k < actual_submitted_tasks; ++k)
    {
        if (k % 10 == 0)
        {
            expected_throwing_tasks++;
        }
        else
        {
            expected_normal_tasks++;
        }
    }

    std::cout << "-------------------- TEST SUMMARY --------------------" << std::endl;
    std::cout << "Tasks initially planned: " << num_tasks << std::endl;
    std::cout << "Tasks actually submitted (via atomic counter): " << actual_submitted_tasks << std::endl;
    std::cout << "Futures added to vector: " << futures.size() << std::endl;
    std::cout << "Futures processed (completed_tasks_count): " << actual_processed_futures << std::endl;
    std::cout << "Counter value (for normal tasks): " << actual_counter_value << std::endl;
    std::cout << "Expected counter value (normal tasks): " << expected_normal_tasks << std::endl;
    std::cout << "Exceptions caught from tasks (via future.get): " << actual_exceptions_from_tasks << std::endl;
    std::cout << "Expected throwing tasks: " << expected_throwing_tasks << std::endl;
    std::cout << "------------------------------------------------------" << std::endl;

    bool pass = true;
    if (actual_counter_value != expected_normal_tasks)
    {
        std::cout << "  [FAIL] Counter mismatch: Expected " << expected_normal_tasks << ", Got " << actual_counter_value << std::endl;
        pass = false;
    }
    // It's possible actual_submitted_tasks < num_tasks if enqueueing itself could fail and be caught,
    // but current ThreadPool::enqueue doesn't throw exceptions that would prevent submission here typically.
    // So, actual_submitted_tasks should ideally be num_tasks.
    if (actual_submitted_tasks != num_tasks)
    {
        std::cout << "  [WARN] Submitted tasks count issue: Planned " << num_tasks << ", Actually Submitted " << actual_submitted_tasks << std::endl;
        // Decide if this is a hard fail depending on ThreadPool::enqueue behavior under stress.
    }
    if (futures.size() != actual_submitted_tasks)
    {
        // This implies an issue with adding futures to the vector or submitted_tasks_atomic_count logic
        std::cout << "  [FAIL] Futures vector size mismatch: Expected " << actual_submitted_tasks
                  << " (based on submitted_tasks_atomic_count), Got " << futures.size() << std::endl;
        pass = false;
    }
    if (actual_processed_futures != futures.size())
    {
        // All futures that were put into the vector should have been processed.
        std::cout << "  [FAIL] Processed futures count mismatch: Expected to process " << futures.size()
                  << " futures, but processed " << actual_processed_futures << std::endl;
        pass = false;
    }
    if (actual_exceptions_from_tasks != expected_throwing_tasks)
    {
        std::cout << "  [FAIL] Exception count mismatch: Expected " << expected_throwing_tasks << ", Got " << actual_exceptions_from_tasks << std::endl;
        pass = false;
    }

    if (pass)
    {
        std::cout << "ThreadPool safety test PASSED!" << std::endl;
        return 0;
    }
    else
    {
        std::cout << "ThreadPool safety test FAILED!" << std::endl;
        return 1;
    }
}