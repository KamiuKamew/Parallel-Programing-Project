#pragma once

#include "../general/utils.h"
#include "../general/op.h"
#include <pthread.h>
#include <vector>
#include <queue>
#include <functional>
#include <unistd.h>

/**
 * @brief 使用NTT优化的多项式乘法
 *
 * @param a 多项式系数
 * @param b 多项式系数
 * @param ab 结果
 * @param n 多项式长度
 * @param p 模数（质数）
 */
template <typename T>
inline void poly_multiply_ntt(T *a, T *b, T *ab, T n, T p, T OMEGA = 3)
{
  // num_threads is now implicitly handled by the global thread pool size.
  // No longer need to pass num_threads to ntt_forward/inverse_mont_pthread.

  using T_mont = T;
  MontMod<T> montMod(p);

  T n_expanded = expand_n(2 * n - 1);
  T *a_expanded = expand_a(a, n, n_expanded);
  T *b_expanded = expand_a(b, n, n_expanded);

  bit_reverse_permute(a_expanded, n_expanded);
  bit_reverse_permute(b_expanded, n_expanded);

  T_mont *a_mont = new T_mont[n_expanded]{};
  T_mont *b_mont = new T_mont[n_expanded]{};
  T_mont *ab_mont = new T_mont[n_expanded]{};
  for (T i = 0; i < n_expanded; ++i)
    a_mont[i] = montMod.from_T(a_expanded[i]);
  for (T i = 0; i < n_expanded; ++i)
    b_mont[i] = montMod.from_T(b_expanded[i]);

  T_mont omega_mont_primitive_root = montMod.from_T(OMEGA);

  ntt_forward_mont_pthread(a_mont, n_expanded, p, omega_mont_primitive_root, montMod);
  ntt_forward_mont_pthread(b_mont, n_expanded, p, omega_mont_primitive_root, montMod);

  for (T i = 0; i < n_expanded; ++i)
    ab_mont[i] = montMod.mul(a_mont[i], b_mont[i]);

  T_mont inv_omega_mont_primitive_root = montMod.inv(omega_mont_primitive_root);
  ntt_inverse_mont_pthread(ab_mont, n_expanded, p, inv_omega_mont_primitive_root, montMod);

  for (T i = 0; i < n_expanded; ++i)
    ab[i] = montMod.to_T(ab_mont[i]);

  delete[] a_expanded;
  delete[] b_expanded;
  delete[] a_mont;
  delete[] b_mont;
  delete[] ab_mont;

  bit_reverse_permute(ab, n_expanded);
}

// Forward declaration for the ThreadPool's static worker function
static void *thread_pool_worker_function(void *context);

class ThreadPool
{
public:
  ThreadPool(size_t num_threads_in_pool) : stop_pool_flag(false), current_active_tasks(0)
  {
    pthread_mutex_init(&task_queue_mutex, nullptr);
    pthread_cond_init(&notify_new_task_condition, nullptr);
    pthread_cond_init(&all_tasks_completed_condition, nullptr);

    worker_pthread_ids.reserve(num_threads_in_pool);
    for (size_t i = 0; i < num_threads_in_pool; ++i)
    {
      pthread_t tid;
      if (pthread_create(&tid, nullptr, thread_pool_worker_function, this) != 0)
      {
        // Handle error: thread creation failed. For simplicity, not throwing.
        // In a real app, this should be handled robustly.
      }
      worker_pthread_ids.push_back(tid);
    }
  }

  ~ThreadPool()
  {
    {
      pthread_mutex_lock(&task_queue_mutex);
      stop_pool_flag = true;
    }
    pthread_cond_broadcast(&notify_new_task_condition); // Wake up all worker threads
    for (pthread_t tid : worker_pthread_ids)
    {
      pthread_join(tid, nullptr);
    }
    pthread_mutex_destroy(&task_queue_mutex);
    pthread_cond_destroy(&notify_new_task_condition);
    pthread_cond_destroy(&all_tasks_completed_condition);
  }

  void enqueue_task(std::function<void()> task_fn)
  {
    pthread_mutex_lock(&task_queue_mutex);
    if (stop_pool_flag)
    {
      pthread_mutex_unlock(&task_queue_mutex);
      // Optionally throw std::runtime_error("enqueue on stopped ThreadPool");
      return;
    }
    internal_task_queue.push(std::move(task_fn));
    current_active_tasks++;
    pthread_mutex_unlock(&task_queue_mutex);
    pthread_cond_signal(&notify_new_task_condition); // Notify one worker thread
  }

  void wait_all_current_tasks_completed()
  {
    pthread_mutex_lock(&task_queue_mutex);
    while (current_active_tasks > 0)
    {
      pthread_cond_wait(&all_tasks_completed_condition, &task_queue_mutex);
    }
    pthread_mutex_unlock(&task_queue_mutex);
  }

  // public for the static C-style worker function to access them via 'this' pointer
public:
  std::vector<pthread_t> worker_pthread_ids;
  std::queue<std::function<void()>> internal_task_queue;
  pthread_mutex_t task_queue_mutex;
  pthread_cond_t notify_new_task_condition;
  pthread_cond_t all_tasks_completed_condition;
  bool stop_pool_flag;
  volatile int current_active_tasks; // Ensure visibility or always protect with mutex
};

// Static worker function for the thread pool
static void *thread_pool_worker_function(void *context)
{
  ThreadPool *pool_instance = static_cast<ThreadPool *>(context);
  while (true)
  {
    std::function<void()> task_to_execute;
    {
      pthread_mutex_lock(&pool_instance->task_queue_mutex);
      while (!pool_instance->stop_pool_flag && pool_instance->internal_task_queue.empty())
      {
        pthread_cond_wait(&pool_instance->notify_new_task_condition, &pool_instance->task_queue_mutex);
      }
      if (pool_instance->stop_pool_flag && pool_instance->internal_task_queue.empty())
      {
        pthread_mutex_unlock(&pool_instance->task_queue_mutex);
        return nullptr; // Exit thread
      }
      task_to_execute = std::move(pool_instance->internal_task_queue.front());
      pool_instance->internal_task_queue.pop();
      pthread_mutex_unlock(&pool_instance->task_queue_mutex);
    }

    task_to_execute(); // Execute the task

    {
      pthread_mutex_lock(&pool_instance->task_queue_mutex);
      pool_instance->current_active_tasks--;
      if (pool_instance->current_active_tasks == 0)
      {
        pthread_cond_signal(&pool_instance->all_tasks_completed_condition);
      }
      pthread_mutex_unlock(&pool_instance->task_queue_mutex);
    }
  }
  return nullptr;
}

// Global constant for number of threads in the pool
static const int NTT_GLOBAL_POOL_THREADS = sysconf(_SC_NPROCESSORS_ONLN) > 0 ? sysconf(_SC_NPROCESSORS_ONLN) : 4;
// Global thread pool instance
static ThreadPool g_ntt_thread_pool(NTT_GLOBAL_POOL_THREADS);

// Structure to pass arguments for a single j-block task
template <typename T>
struct NttPthreadTaskArgs
{
  T *a_mont_ptr;
  T n_length;
  T current_mid_val;
  T Wn_mont_val_for_mid;
  MontMod<T> *mont_mod_obj_ptr;
  T j_block_idx; // The specific j index for this task
};

// Helper function to execute the core logic for a single j-block
template <typename T>
void execute_ntt_j_block_computation(NttPthreadTaskArgs<T> args)
{
  T w_mont_iter = args.mont_mod_obj_ptr->from_u32(1);
  // The k-loop processes elements for the given j_block_idx
  for (T k = 0; k < args.current_mid_val; ++k)
  {
    // Ensure indices are within bounds, especially j_block_idx + k + current_mid_val
    if (args.j_block_idx + k + args.current_mid_val >= args.n_length)
    {
      // This check might be optimistic if n_length is not a power of 2.
      // NTT typically assumes power-of-2 lengths for full blocks.
      // If partial blocks are possible and need specific handling:
      // T idx1 = args.j_block_idx + k;
      // T idx2 = args.j_block_idx + k + args.current_mid_val;
      // if (idx2 >= args.n_length) continue; // or break, depending on logic for partial blocks
      continue;
    }
    T x_mont = args.a_mont_ptr[args.j_block_idx + k];
    T y_mont = args.mont_mod_obj_ptr->mul(w_mont_iter, args.a_mont_ptr[args.j_block_idx + k + args.current_mid_val]);
    args.a_mont_ptr[args.j_block_idx + k] = args.mont_mod_obj_ptr->add(x_mont, y_mont);
    args.a_mont_ptr[args.j_block_idx + k + args.current_mid_val] = args.mont_mod_obj_ptr->sub(x_mont, y_mont);
    w_mont_iter = args.mont_mod_obj_ptr->mul(w_mont_iter, args.Wn_mont_val_for_mid);
  }
}

// Pthread version of ntt_forward_mont using the global thread pool
template <typename T>
void ntt_forward_mont_pthread(T *a_mont_arr, T n_len, T p_mod, T omega_mont_root, MontMod<T> &montMod)
{
  for (T mid = 1; mid < n_len; mid <<= 1)
  {
    T Wn_val_current_mid = montMod.pow(omega_mont_root, (p_mod - 1) / (mid << 1));
    T j_loop_increment = (mid << 1);
    T total_j_blocks_for_mid = (n_len / j_loop_increment);

    if (total_j_blocks_for_mid == 0 && n_len > 0 && mid < n_len)
    {
      total_j_blocks_for_mid = 1;
    }

    // If not enough blocks to parallelize effectively, or if n_len is small.
    // The threshold (e.g., < 2 blocks) can be tuned.
    if (total_j_blocks_for_mid < 2 || n_len < (mid << 2) /* Heuristic for small N */)
    { // Serial execution
      for (T j = 0; j < n_len; j += j_loop_increment)
      {
        NttPthreadTaskArgs<T> task_args_serial = {
            a_mont_arr, n_len, mid, Wn_val_current_mid, &montMod, j};
        execute_ntt_j_block_computation(task_args_serial);
      }
    }
    else
    { // Parallel execution using the global thread pool
      for (T j_block_num = 0; j_block_num < total_j_blocks_for_mid; ++j_block_num)
      {
        T current_j_for_task = j_block_num * j_loop_increment;
        NttPthreadTaskArgs<T> task_args_parallel = {
            a_mont_arr, n_len, mid, Wn_val_current_mid, &montMod, current_j_for_task};
        // Capture task_args_parallel by value for the lambda
        g_ntt_thread_pool.enqueue_task([task_args_parallel]()
                                       { execute_ntt_j_block_computation(task_args_parallel); });
      }
      g_ntt_thread_pool.wait_all_current_tasks_completed();
    }
  }
}

// Pthread version of ntt_inverse_mont using the global thread pool
template <typename T>
void ntt_inverse_mont_pthread(T *a_mont_arr, T n_len, T p_mod, T inv_omega_mont_root, MontMod<T> &montMod)
{
  for (T mid = 1; mid < n_len; mid <<= 1)
  {
    T Wn_val_current_mid = montMod.pow(inv_omega_mont_root, (p_mod - 1) / (mid << 1));
    T j_loop_increment = (mid << 1);
    T total_j_blocks_for_mid = (n_len / j_loop_increment);

    if (total_j_blocks_for_mid == 0 && n_len > 0 && mid < n_len)
    {
      total_j_blocks_for_mid = 1;
    }

    if (total_j_blocks_for_mid < 2 || n_len < (mid << 2))
    { // Serial execution
      for (T j = 0; j < n_len; j += j_loop_increment)
      {
        NttPthreadTaskArgs<T> task_args_serial = {
            a_mont_arr, n_len, mid, Wn_val_current_mid, &montMod, j};
        execute_ntt_j_block_computation(task_args_serial);
      }
    }
    else
    { // Parallel execution using the global thread pool
      for (T j_block_num = 0; j_block_num < total_j_blocks_for_mid; ++j_block_num)
      {
        T current_j_for_task = j_block_num * j_loop_increment;
        NttPthreadTaskArgs<T> task_args_parallel = {
            a_mont_arr, n_len, mid, Wn_val_current_mid, &montMod, current_j_for_task};
        g_ntt_thread_pool.enqueue_task([task_args_parallel]()
                                       { execute_ntt_j_block_computation(task_args_parallel); });
      }
      g_ntt_thread_pool.wait_all_current_tasks_completed();
    }
  }

  // Final scaling for inverse NTT
  if (n_len > 0)
  { // Avoid division by zero or issues if n_len is 0
    T n_len_mont = montMod.from_u32(n_len);
    T n_inv_mont = montMod.inv(n_len_mont);
    for (T i = 0; i < n_len; ++i)
    {
      a_mont_arr[i] = montMod.mul(a_mont_arr[i], n_inv_mont);
    }
  }
}