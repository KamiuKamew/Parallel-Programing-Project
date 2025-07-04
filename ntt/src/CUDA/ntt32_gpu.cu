#include <algorithm>
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "ntt32_gpu.h"
#include <cuda_runtime.h>

#define CUDA_CHECK(err)                                                        \
  {                                                                            \
    cudaError_t err_ = (err);                                                  \
    if (err_ != cudaSuccess) {                                                 \
      std::cerr << "CUDA error in " << __FILE__ << " at line " << __LINE__     \
                << ": " << cudaGetErrorString(err_) << std::endl;              \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  }

// 基础NTT kernel - 修复线程索引计算
__global__ void ntt_kernel_basic(u32 *a, const u32 *twiddles, int len, int n,
                                 u32 mod, bool invert) {
  int tidx = blockIdx.x * blockDim.x + threadIdx.x;

  // 计算当前线程负责的蝶形变换
  int butterflies_per_stage = n / len * (len / 2);
  if (tidx >= butterflies_per_stage)
    return;

  int group_size = len;
  int butterflies_per_group = len / 2;
  int group_idx = tidx / butterflies_per_group;
  int butterfly_idx = tidx % butterflies_per_group;

  int i = group_idx * group_size + butterfly_idx;
  int j = i + butterflies_per_group;

  if (i < n && j < n) {
    u32 w = twiddles[butterfly_idx];
    u32 u = a[i];
    u32 v = mul_mod(a[j], w, mod);

    a[i] = add_mod(u, v, mod);
    a[j] = sub_mod(u, v, mod);
  }
}

// Barrett优化版本
template <typename Reducer>
__global__ void ntt_kernel_optimized(u32 *a, const u32 *twiddles, int len,
                                     int n, const Reducer reducer,
                                     bool invert) {
  int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  int butterfly_grp_idx = tidx / (len / 2);
  int butterfly_idx_in_grp = tidx % (len / 2);
  int i = butterfly_grp_idx * len + butterfly_idx_in_grp;

  if (i < n && i + len / 2 < n) {
    u32 w = twiddles[butterfly_idx_in_grp];
    u32 u = a[i];
    u32 v = reducer.multiply(a[i + len / 2], w);

    u32 sum = u + v;
    a[i] = (sum >= reducer.mod) ? (sum - reducer.mod) : sum;

    u32 diff = (u >= v) ? (u - v) : (u + reducer.mod - v);
    a[i + len / 2] = diff;
  }
}

// 点乘kernel
__global__ void pointwise_mult_kernel(u32 *out, const u32 *in1, const u32 *in2,
                                      int n, u32 mod) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    out[idx] = mul_mod(in1[idx], in2[idx], mod);
  }
}

template <typename Reducer>
__global__ void pointwise_mult_kernel_optimized(u32 *out, const u32 *in1,
                                                const u32 *in2, int n,
                                                const Reducer reducer) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    out[idx] = reducer.multiply(in1[idx], in2[idx]);
  }
}

// 最终缩放kernel
__global__ void final_scaling_kernel(u32 *a, int n, u32 n_inv, u32 mod) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    a[idx] = mul_mod(a[idx], n_inv, mod);
  }
}

template <typename Reducer>
__global__ void final_scaling_kernel_optimized(u32 *a, int n, u32 n_inv,
                                               const Reducer reducer) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    a[idx] = reducer.multiply(a[idx], n_inv);
  }
}

// Montgomery相关kernel
template <typename Reducer>
__global__ void to_mont_kernel(u32 *a, int n, const Reducer reducer) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    a[idx] = reducer.to_mont(a[idx]);
  }
}

template <typename Reducer>
__global__ void from_mont_kernel(u32 *a, int n, const Reducer reducer) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    a[idx] = reducer.from_mont(a[idx]);
  }
}

// 添加bit-reverse kernel
__global__ void bit_reverse_kernel(u32 *a, int n) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= n)
    return;

  int lg_n = 0;
  int temp_n = n;
  while (temp_n > 1) {
    lg_n++;
    temp_n >>= 1;
  }

  int j = 0;
  for (int k = 0; k < lg_n; ++k) {
    if (idx & (1 << k)) {
      j |= (1 << (lg_n - 1 - k));
    }
  }

  if (idx < j) {
    u32 temp = a[idx];
    a[idx] = a[j];
    a[j] = temp;
  }
}

// 单线程GPU NTT kernel（用于调试）
__global__ void ntt_kernel_single_thread(u32 *a, int n, u32 mod, u32 root,
                                         bool invert) {
  if (threadIdx.x != 0 || blockIdx.x != 0)
    return;

  // Bit-reverse
  for (int i = 1, j = 0; i < n; i++) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) {
      u32 temp = a[i];
      a[i] = a[j];
      a[j] = temp;
    }
  }

  // NTT
  for (int len = 2; len <= n; len <<= 1) {
    u32 wlen = pow_mod(root, (mod - 1) / len, mod);
    if (invert)
      wlen = inv_mod(wlen, mod);

    for (int i = 0; i < n; i += len) {
      u32 w = 1;
      for (int j = 0; j < len / 2; j++) {
        u32 u = a[i + j];
        u32 v = mul_mod(a[i + j + len / 2], w, mod);
        a[i + j] = add_mod(u, v, mod);
        a[i + j + len / 2] = sub_mod(u, v, mod);
        w = mul_mod(w, wlen, mod);
      }
    }
  }

  if (invert) {
    u32 n_inv = inv_mod(n, mod);
    for (int i = 0; i < n; i++) {
      a[i] = mul_mod(a[i], n_inv, mod);
    }
  }
}

// 单线程GPU NTT实现
void ntt_gpu_single_thread(u32 *d_a, int n, bool invert, u32 mod, u32 root) {
  ntt_kernel_single_thread<<<1, 1>>>(d_a, n, mod, root, invert);
  CUDA_CHECK(cudaGetLastError());
}

// 优化的多线程NTT kernel
__global__ void ntt_kernel_optimized_mt(u32 *a, int n, u32 mod, u32 root,
                                        bool invert) {
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  // 每个线程处理一个数组元素的bit-reverse
  if (tid < n) {
    int lg_n = 0;
    int temp_n = n;
    while (temp_n > 1) {
      lg_n++;
      temp_n >>= 1;
    }

    int j = 0;
    for (int k = 0; k < lg_n; ++k) {
      if (tid & (1 << k)) {
        j |= (1 << (lg_n - 1 - k));
      }
    }

    // 只有当 tid < j 时才交换，避免重复交换
    if (tid < j) {
      u32 temp = a[tid];
      a[tid] = a[j];
      a[j] = temp;
    }
  }

  __syncthreads();

  // NTT变换
  for (int len = 2; len <= n; len <<= 1) {
    u32 wlen = pow_mod(root, (mod - 1) / len, mod);
    if (invert)
      wlen = inv_mod(wlen, mod);

    int num_groups = n / len;
    int group_id = tid / (len / 2);
    int elem_id = tid % (len / 2);

    if (group_id < num_groups) {
      int base = group_id * len;
      int i = base + elem_id;
      int j = i + len / 2;

      u32 w = pow_mod(wlen, elem_id, mod);
      u32 u = a[i];
      u32 v = mul_mod(a[j], w, mod);

      a[i] = add_mod(u, v, mod);
      a[j] = sub_mod(u, v, mod);
    }

    __syncthreads();
  }

  // 逆变换时的缩放
  if (invert && tid < n) {
    u32 n_inv = inv_mod(n, mod);
    a[tid] = mul_mod(a[tid], n_inv, mod);
  }
}

// 多线程GPU NTT实现
void ntt_gpu_optimized_mt(u32 *d_a, int n, bool invert, u32 mod, u32 root) {
  int threads = min(n, 1024); // 使用单个block，最多1024线程
  int blocks = 1;
  ntt_kernel_optimized_mt<<<blocks, threads>>>(d_a, n, mod, root, invert);
  CUDA_CHECK(cudaGetLastError());
}

// GPU NTT实现 - 基础版本
void ntt_gpu_basic(u32 *d_a, int n, bool invert, u32 mod, u32 root) {
  int threads_per_block = 256;
  int blocks = (n + threads_per_block - 1) / threads_per_block;

  // 前向NTT需要先bit-reverse
  if (!invert) {
    bit_reverse_kernel<<<blocks, threads_per_block>>>(d_a, n);
    CUDA_CHECK(cudaGetLastError());
  }

  for (int len = 2; len <= n; len <<= 1) {
    std::vector<u32> h_twiddles(len / 2);
    u32 wlen = pow_mod(root, (mod - 1) / len, mod); // w_len = g^((p-1)/len)
    if (invert)
      wlen = inv_mod(wlen, mod);

    h_twiddles[0] = 1;
    for (int j = 1; j < len / 2; j++) {
      h_twiddles[j] = mul_mod(h_twiddles[j - 1], wlen, mod);
    }

    u32 *d_twiddles;
    CUDA_CHECK(cudaMalloc((void **)&d_twiddles, (len / 2) * sizeof(u32)));
    CUDA_CHECK(cudaMemcpy(d_twiddles, h_twiddles.data(),
                          (len / 2) * sizeof(u32), cudaMemcpyHostToDevice));

    int num_threads = n / len * (len / 2); // 每个长度len的组有len/2个蝶形
    int num_blocks = (num_threads + threads_per_block - 1) / threads_per_block;
    ntt_kernel_basic<<<num_blocks, threads_per_block>>>(d_a, d_twiddles, len, n,
                                                        mod, invert);
    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaFree(d_twiddles));
  }

  // 逆向NTT后需要bit-reverse恢复自然顺序
  if (invert) {
    bit_reverse_kernel<<<blocks, threads_per_block>>>(d_a, n);
    CUDA_CHECK(cudaGetLastError());
  }
}

// GPU NTT实现 - Barrett优化版本
template <typename Reducer>
void ntt_gpu_optimized(u32 *d_a, int n, bool invert, const Reducer &reducer) {
  int threads_per_block = 256;

  // 预计算所有旋转因子
  std::vector<u32> h_all_twiddles;
  size_t total_twiddles = 0;
  for (int len = 2; len <= n; len <<= 1) {
    total_twiddles += len / 2;
  }
  h_all_twiddles.reserve(total_twiddles);

  u32 root = pow_mod(3, (reducer.mod - 1) / n, reducer.mod);
  if (invert) {
    root = inv_mod(root, reducer.mod);
  }

  for (int len = 2; len <= n; len <<= 1) {
    u32 wlen_base = pow_mod(root, n / len, reducer.mod);
    u32 w = 1;
    for (int j = 0; j < len / 2; j++) {
      h_all_twiddles.push_back(w);
      w = mul_mod(w, wlen_base, reducer.mod);
    }
  }

  u32 *d_all_twiddles;
  CUDA_CHECK(
      cudaMalloc((void **)&d_all_twiddles, total_twiddles * sizeof(u32)));
  CUDA_CHECK(cudaMemcpy(d_all_twiddles, h_all_twiddles.data(),
                        total_twiddles * sizeof(u32), cudaMemcpyHostToDevice));

  size_t twiddle_offset = 0;
  for (int len = 2; len <= n; len <<= 1) {
    u32 *d_twiddles_stage = d_all_twiddles + twiddle_offset;

    int num_threads = n / 2;
    int num_blocks = (num_threads + threads_per_block - 1) / threads_per_block;

    ntt_kernel_optimized<Reducer><<<num_blocks, threads_per_block>>>(
        d_a, d_twiddles_stage, len, n, reducer, invert);
    CUDA_CHECK(cudaGetLastError());

    twiddle_offset += len / 2;
  }

  CUDA_CHECK(cudaFree(d_all_twiddles));
}

// 添加 BasicReducer32（仅使用 mul_mod）
struct BasicReducer32 {
  u32 mod;
  __host__ __device__ BasicReducer32(u32 m = 0) : mod(m) {}
  __host__ __device__ u32 multiply(u32 a, u32 b) const {
    return mul_mod(a, b, mod);
  }
};

// 共享内存版单 stage kernel，适用于 len <= 2048
__global__ void ntt_stage_shared(u32 *a, const u32 *twiddles, int len, int n,
                                 u32 mod) {
  extern __shared__ u32 s_twiddles[]; // len/2 元素
  int tid = threadIdx.x;
  int half = len >> 1;
  if (tid < half) {
    s_twiddles[tid] = twiddles[tid];
  }
  __syncthreads();

  int group = blockIdx.x; // 当前 block 对应的蝶形组编号
  int idx = tid;
  while (idx < half) {
    int i = group * len + idx;
    int j = i + half;
    if (j < n) {
      u32 w = s_twiddles[idx];
      u32 u = a[i];
      u32 v = mul_mod(a[j], w, mod);
      a[i] = add_mod(u, v, mod);
      a[j] = sub_mod(u, v, mod);
    }
    idx += blockDim.x; // 如果 blockDim.x 小于 half，让线程继续处理剩余蝶形
  }
}

// 使用共享内存 + 多 block 的 NTT（仅 forward/inverse 变换，不做缩放）
void ntt_gpu_shared(u32 *d_a, int n, bool invert, u32 mod, u32 root) {
  constexpr int SMEM_THRESHOLD = 2048; // 只有 len<=2048 时才使用共享内存 kernel
  // Bit-reverse copy on GPU
  int threads_per_block = 256;
  int blocks_rev = (n + threads_per_block - 1) / threads_per_block;
  if (!invert) {
    bit_reverse_kernel<<<blocks_rev, threads_per_block>>>(d_a, n);
    CUDA_CHECK(cudaGetLastError());
  }

  // 预先在 host 端一次性计算所有 twiddle 并拷贝上去（与 ntt_gpu_optimized
  // 相同）
  std::vector<u32> h_all_twiddles;
  size_t total_twiddles = 0;
  for (int len = 2; len <= n; len <<= 1)
    total_twiddles += len / 2;
  h_all_twiddles.reserve(total_twiddles);

  u32 omega = pow_mod(root, (mod - 1) / n, mod);
  if (invert)
    omega = inv_mod(omega, mod);

  for (int len = 2; len <= n; len <<= 1) {
    u32 step = n / len;
    u32 w = 1;
    for (int j = 0; j < len / 2; ++j) {
      h_all_twiddles.push_back(w);
      w = mul_mod(w, pow_mod(omega, step, mod), mod);
    }
  }
  u32 *d_all_twiddles;
  CUDA_CHECK(
      cudaMalloc((void **)&d_all_twiddles, total_twiddles * sizeof(u32)));
  CUDA_CHECK(cudaMemcpy(d_all_twiddles, h_all_twiddles.data(),
                        total_twiddles * sizeof(u32), cudaMemcpyHostToDevice));

  size_t offset = 0;
  for (int len = 2; len <= n; len <<= 1) {
    u32 *d_twiddles_stage = d_all_twiddles + offset;
    int half = len >> 1;
    int blocks = n / len; // 每个 block 负责一个蝶形组
    if (len <= SMEM_THRESHOLD) {
      int threads = (half < 1024) ? half : 1024;
      size_t shared_bytes = half * sizeof(u32);
      ntt_stage_shared<<<blocks, threads, shared_bytes>>>(d_a, d_twiddles_stage,
                                                          len, n, mod);
    } else {
      // 使用已有的全局内存 kernel（与 ntt_kernel_basic 相同逻辑）
      int total_threads = n / 2;
      int threads = 256;
      int blocks_glb = (total_threads + threads - 1) / threads;
      ntt_kernel_basic<<<blocks_glb, threads>>>(d_a, d_twiddles_stage, len, n,
                                                mod, invert);
    }
    CUDA_CHECK(cudaGetLastError());
    offset += half;
  }

  if (invert) {
    bit_reverse_kernel<<<blocks_rev, threads_per_block>>>(d_a, n);
    CUDA_CHECK(cudaGetLastError());
  }

  CUDA_CHECK(cudaFree(d_all_twiddles));
}

// 主要接口函数
std::vector<u32> multiply_ntt_gpu32(std::vector<u32> &poly1,
                                    std::vector<u32> &poly2, u32 mod,
                                    u32 primitive_root,
                                    const std::string &method) {

  int n1 = poly1.size();
  int n2 = poly2.size();
  if (n1 == 0 || n2 == 0)
    return {};

  int target_len = n1 + n2 - 1;
  int n = 1;
  while (n < target_len)
    n <<= 1;

  poly1.resize(n);
  poly2.resize(n);

  u32 *d_p1, *d_p2;
  CUDA_CHECK(cudaMalloc((void **)&d_p1, n * sizeof(u32)));
  CUDA_CHECK(cudaMalloc((void **)&d_p2, n * sizeof(u32)));
  CUDA_CHECK(
      cudaMemcpy(d_p1, poly1.data(), n * sizeof(u32), cudaMemcpyHostToDevice));
  CUDA_CHECK(
      cudaMemcpy(d_p2, poly2.data(), n * sizeof(u32), cudaMemcpyHostToDevice));

  int threads = 256;
  int blocks = (n + threads - 1) / threads;

  if (method == "barrett") {
    BarrettReducer32 br(mod);
    ntt_gpu_optimized<BarrettReducer32>(d_p1, n, false, br);
    ntt_gpu_optimized<BarrettReducer32>(d_p2, n, false, br);
    pointwise_mult_kernel_optimized<BarrettReducer32>
        <<<blocks, threads>>>(d_p1, d_p1, d_p2, n, br);
    ntt_gpu_optimized<BarrettReducer32>(d_p1, n, true, br);
    u32 n_inv = inv_mod(n, mod);
    final_scaling_kernel_optimized<BarrettReducer32>
        <<<blocks, threads>>>(d_p1, n, n_inv, br);
  } else if (method == "montgomery") {
    MontgomeryReducer32 mr(mod);

    to_mont_kernel<MontgomeryReducer32><<<blocks, threads>>>(d_p1, n, mr);
    to_mont_kernel<MontgomeryReducer32><<<blocks, threads>>>(d_p2, n, mr);

    ntt_gpu_optimized<MontgomeryReducer32>(d_p1, n, false, mr);
    ntt_gpu_optimized<MontgomeryReducer32>(d_p2, n, false, mr);

    pointwise_mult_kernel_optimized<MontgomeryReducer32>
        <<<blocks, threads>>>(d_p1, d_p1, d_p2, n, mr);

    ntt_gpu_optimized<MontgomeryReducer32>(d_p1, n, true, mr);

    u32 n_inv = inv_mod(n, mod);
    u32 n_inv_mont = mr.to_mont(n_inv);
    final_scaling_kernel_optimized<MontgomeryReducer32>
        <<<blocks, threads>>>(d_p1, n, n_inv_mont, mr);

    from_mont_kernel<MontgomeryReducer32><<<blocks, threads>>>(d_p1, n, mr);

  } else if (method == "fast") {
    ntt_gpu_shared(d_p1, n, false, mod, primitive_root);
    ntt_gpu_shared(d_p2, n, false, mod, primitive_root);
    pointwise_mult_kernel<<<blocks, threads>>>(d_p1, d_p1, d_p2, n, mod);
    ntt_gpu_shared(d_p1, n, true, mod, primitive_root);
    u32 n_inv = inv_mod(n, mod);
    final_scaling_kernel<<<blocks, threads>>>(d_p1, n, n_inv, mod);
  } else { // basic
    if (method == "optimized") {
      ntt_gpu_optimized_mt(d_p1, n, false, mod, primitive_root);
      ntt_gpu_optimized_mt(d_p2, n, false, mod, primitive_root);
    } else {
      ntt_gpu_single_thread(d_p1, n, false, mod, primitive_root);
      ntt_gpu_single_thread(d_p2, n, false, mod, primitive_root);
    }
    pointwise_mult_kernel<<<blocks, threads>>>(d_p1, d_p1, d_p2, n, mod);
    if (method == "optimized") {
      ntt_gpu_optimized_mt(d_p1, n, true, mod, primitive_root);
    } else {
      ntt_gpu_single_thread(d_p1, n, true, mod, primitive_root);
    }
    // 多线程和单线程版本都已经包含缩放，不需要额外调用
  }

  std::vector<u32> result(n);
  CUDA_CHECK(
      cudaMemcpy(result.data(), d_p1, n * sizeof(u32), cudaMemcpyDeviceToHost));

  CUDA_CHECK(cudaFree(d_p1));
  CUDA_CHECK(cudaFree(d_p2));

  result.resize(target_len);
  return result;
}

// 调试用：简单的CPU版本NTT（用于对比）
std::vector<u32> multiply_ntt_cpu32_debug(std::vector<u32> &poly1,
                                          std::vector<u32> &poly2, u32 mod,
                                          u32 primitive_root) {
  int n1 = poly1.size();
  int n2 = poly2.size();
  if (n1 == 0 || n2 == 0)
    return {};

  int target_len = n1 + n2 - 1;
  int n = 1;
  while (n < target_len)
    n <<= 1;

  std::vector<u32> a(poly1.begin(), poly1.end());
  std::vector<u32> b(poly2.begin(), poly2.end());
  a.resize(n);
  b.resize(n);

  // 简单的CPU NTT
  auto ntt_cpu = [&](std::vector<u32> &vec, bool invert) {
    int n = vec.size();

    // Bit-reverse
    for (int i = 1, j = 0; i < n; i++) {
      int bit = n >> 1;
      for (; j & bit; bit >>= 1) {
        j ^= bit;
      }
      j ^= bit;
      if (i < j)
        std::swap(vec[i], vec[j]);
    }

    // NTT
    for (int len = 2; len <= n; len <<= 1) {
      u32 wlen = pow_mod(primitive_root, (mod - 1) / len, mod);
      if (invert)
        wlen = inv_mod(wlen, mod);

      for (int i = 0; i < n; i += len) {
        u32 w = 1;
        for (int j = 0; j < len / 2; j++) {
          u32 u = vec[i + j];
          u32 v = mul_mod(vec[i + j + len / 2], w, mod);
          vec[i + j] = add_mod(u, v, mod);
          vec[i + j + len / 2] = sub_mod(u, v, mod);
          w = mul_mod(w, wlen, mod);
        }
      }
    }

    if (invert) {
      u32 n_inv = inv_mod(n, mod);
      for (auto &x : vec) {
        x = mul_mod(x, n_inv, mod);
      }
    }
  };

  ntt_cpu(a, false);
  ntt_cpu(b, false);

  for (int i = 0; i < n; i++) {
    a[i] = mul_mod(a[i], b[i], mod);
  }

  ntt_cpu(a, true);

  a.resize(target_len);
  return a;
}