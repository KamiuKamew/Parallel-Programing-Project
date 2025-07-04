#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <cuda_runtime.h>

#ifdef __NVCC__
#define CUDA_DEVICE_CODE
#endif

#ifdef __GNUC__
using int128 = __int128;
#else
using int128 = long long;
#endif

using ll = long long;

__host__ __device__ inline ll power(ll base, ll exp, ll mod) {
  ll res = 1;
  base %= mod;
  while (exp > 0) {
    if (exp % 2 == 1)
      res = (int128)res * base % mod;
    base = (int128)base * base % mod;
    exp /= 2;
  }
  return res;
}

__host__ __device__ inline ll modInverse(ll n, ll mod) {
  return power(n, mod - 2, mod);
}

__host__ __device__ inline ll extended_gcd_device(ll a, ll b, ll &x, ll &y) {
  if (a == 0) {
    x = 0;
    y = 1;
    return b;
  }
  ll x1, y1;
  ll gcd = extended_gcd_device(b % a, a, x1, y1);
  x = y1 - (b / a) * x1;
  y = x1;
  return gcd;
}

struct BarrettReducer {
  ll mod;
  unsigned __int128 mu;

  __host__ __device__ BarrettReducer(ll m) : mod(m), mu(0) {
    if (m != 0) {
      mu = (unsigned __int128)-1 / m + 1;
    }
  }

  __host__ __device__ ll multiply(ll a, ll b) const {
    if (mod == 0)
      return (ll)((unsigned __int128)a * b);
    unsigned __int128 prod = (unsigned __int128)a * b;
    unsigned __int128 q = (prod * mu) >> 64;
    ll r = (ll)(prod - q * mod);
    if (r < 0)
      r += mod;
    if (r >= mod)
      r -= mod;
    return r;
  }
};

struct MontgomeryReducer {
  ll mod;
  ll r_inv;
  ll mod_prime;
  static constexpr int R_BITS = 62;
  static constexpr ll R = 1LL << R_BITS;

  __host__ __device__ MontgomeryReducer(ll m = 0)
      : mod(m), r_inv(0), mod_prime(0) {
    if (mod == 0)
      return;

    ll y;
    extended_gcd_device(R, mod, r_inv, y);
    r_inv = (r_inv % mod + mod) % mod;

    mod_prime = 0;
    ll t_inv = 0;
    for (int i = 0; i < R_BITS; ++i) {
      if (!(t_inv & 1)) {
        t_inv += mod;
        mod_prime |= (1LL << i);
      }
      t_inv >>= 1;
    }
  }

  __host__ __device__ ll to_mont(ll a) const {
    return (int128(a) << R_BITS) % mod;
  }

  __host__ __device__ ll from_mont(ll a_mont) const {
    ll m = ((int128(a_mont) * mod_prime) & (R - 1));
    ll t = (a_mont + m * mod) >> R_BITS;
    if (t >= mod)
      return t - mod;
    return t;
  }

  __host__ __device__ ll multiply(ll a_mont, ll b_mont) const {
    int128 t = (int128)a_mont * b_mont;
    ll m = ((t & (R - 1)) * mod_prime) & (R - 1);
    ll u = (t + m * mod) >> R_BITS;
    if (u >= mod)
      return u - mod;
    return u;
  }
};

#define CUDA_CHECK(err)                                                        \
  {                                                                            \
    cudaError_t err_ = (err);                                                  \
    if (err_ != cudaSuccess) {                                                 \
      std::cerr << "CUDA error in " << __FILE__ << " at line " << __LINE__     \
                << ": " << cudaGetErrorString(err_) << std::endl;              \
      exit(EXIT_FAILURE);                                                      \
    }                                                                          \
  }

__global__ void pointwise_mult_kernel(ll *out, const ll *in1, const ll *in2,
                                      int n, ll mod) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    out[idx] = ((int128)in1[idx] * in2[idx]) % mod;
  }
}

__global__ void final_scaling_kernel(ll *a, int n, ll n_inv, ll mod) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    a[idx] = ((int128)a[idx] * n_inv) % mod;
  }
}

__global__ void ntt_kernel_naive(ll *a, const ll *twiddles, int len, int n,
                                 ll mod) {
  int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  int butterfly_grp_idx = tidx / (len / 2);
  int butterfly_idx_in_grp = tidx % (len / 2);
  int i = butterfly_grp_idx * len + butterfly_idx_in_grp;

  if (i < n) {
    ll w = twiddles[butterfly_idx_in_grp];
    ll u = a[i];
    ll v = ((int128)a[i + len / 2] * w) % mod;
    a[i] = (u + v) % mod;
    a[i + len / 2] = (u - v + mod) % mod;
  }
}

template <typename Reducer>
__global__ void pointwise_mult_kernel_optimized(ll *out, const ll *in1,
                                                const ll *in2, int n,
                                                const Reducer reducer) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    out[idx] = reducer.multiply(in1[idx], in2[idx]);
  }
}

template <typename Reducer>
__global__ void final_scaling_kernel_optimized(ll *a, int n, ll n_inv,
                                               const Reducer reducer) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    a[idx] = reducer.multiply(a[idx], n_inv);
  }
}

template <typename Reducer>
__global__ void ntt_kernel_optimized(ll *a, const ll *twiddles, int len, int n,
                                     const Reducer reducer) {
  int tidx = blockIdx.x * blockDim.x + threadIdx.x;
  int butterfly_grp_idx = tidx / (len / 2);
  int butterfly_idx_in_grp = tidx % (len / 2);
  int i = butterfly_grp_idx * len + butterfly_idx_in_grp;

  if (i < n) {
    ll w = twiddles[butterfly_idx_in_grp];
    ll u = a[i];
    ll v = reducer.multiply(a[i + len / 2], w);

    ll sum = u + v;
    a[i] = (sum >= reducer.mod) ? (sum - reducer.mod) : sum;

    ll diff = u - v;
    a[i + len / 2] = (diff < 0) ? (diff + reducer.mod) : diff;
  }
}

template <typename Reducer>
__global__ void from_mont_kernel(ll *a, int n, const Reducer reducer) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    a[idx] = reducer.from_mont(a[idx]);
  }
}

template <typename Reducer>
__global__ void to_mont_kernel(ll *a, int n, const Reducer reducer) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < n) {
    a[idx] = reducer.to_mont(a[idx]);
  }
}

void bit_reverse_copy(std::vector<ll> &a) {
  int n = a.size();
  for (int i = 1, j = 0; i < n; i++) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) {
      j ^= bit;
    }
    j ^= bit;
    if (i < j) {
      std::swap(a[i], a[j]);
    }
  }
}

void ntt_gpu_naive(ll *d_a, int n, bool invert, ll mod, ll root) {
  int threads_per_block = 256;
  for (int len = 2; len <= n; len <<= 1) {
    std::vector<ll> h_twiddles(len / 2);
    ll wlen = power(root, (mod - 1) / len, mod);
    if (invert)
      wlen = modInverse(wlen, mod);

    h_twiddles[0] = 1;
    for (int j = 1; j < len / 2; j++) {
      h_twiddles[j] = ((int128)h_twiddles[j - 1] * wlen) % mod;
    }

    ll *d_twiddles;
    CUDA_CHECK(cudaMalloc((void **)&d_twiddles, (len / 2) * sizeof(ll)));
    CUDA_CHECK(cudaMemcpy(d_twiddles, h_twiddles.data(), (len / 2) * sizeof(ll),
                          cudaMemcpyHostToDevice));

    int num_threads = n / 2;
    int num_blocks = (num_threads + threads_per_block - 1) / threads_per_block;
    ntt_kernel_naive<<<num_blocks, threads_per_block>>>(d_a, d_twiddles, len, n,
                                                        mod);
    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaFree(d_twiddles));
  }
}

template <typename Reducer>
void ntt_gpu_optimized(ll *d_a, int n, bool invert, const Reducer &reducer) {
  int threads_per_block = 256;

  std::vector<ll> h_all_twiddles;
  size_t total_twiddles = 0;
  for (int len = 2; len <= n; len <<= 1) {
    total_twiddles += len / 2;
  }
  h_all_twiddles.reserve(total_twiddles);

  ll root = power(3, (reducer.mod - 1) / n, reducer.mod);
  if (invert) {
    root = modInverse(root, reducer.mod);
  }

  for (int len = 2; len <= n; len <<= 1) {
    ll wlen_base = power(root, n / len, reducer.mod);
    ll w = 1;
    for (int j = 0; j < len / 2; j++) {
      h_all_twiddles.push_back(w);
      w = ((int128)w * wlen_base) % reducer.mod;
    }
  }

  ll *d_all_twiddles;
  CUDA_CHECK(cudaMalloc((void **)&d_all_twiddles, total_twiddles * sizeof(ll)));
  CUDA_CHECK(cudaMemcpy(d_all_twiddles, h_all_twiddles.data(),
                        total_twiddles * sizeof(ll), cudaMemcpyHostToDevice));

  size_t twiddle_offset = 0;
  for (int len = 2; len <= n; len <<= 1) {
    ll *d_twiddles_stage = d_all_twiddles + twiddle_offset;

    int num_threads = n / 2;
    int num_blocks = (num_threads + threads_per_block - 1) / threads_per_block;

    ntt_kernel_optimized<Reducer><<<num_blocks, threads_per_block>>>(
        d_a, d_twiddles_stage, len, n, reducer);
    CUDA_CHECK(cudaGetLastError());

    twiddle_offset += len / 2;
  }

  CUDA_CHECK(cudaFree(d_all_twiddles));
}

template <typename Reducer>
void ntt_gpu_montgomery(ll *d_a, int n, bool invert, const Reducer &reducer) {
  int threads_per_block = 256;

  std::vector<ll> h_all_twiddles;
  size_t total_twiddles = 0;
  for (int len = 2; len <= n; len <<= 1) {
    total_twiddles += len / 2;
  }
  h_all_twiddles.reserve(total_twiddles);

  ll root = power(3, (reducer.mod - 1) / n, reducer.mod);
  if (invert) {
    root = modInverse(root, reducer.mod);
  }

  for (int len = 2; len <= n; len <<= 1) {
    ll wlen_base = power(root, n / len, reducer.mod);
    ll w = 1;
    for (int j = 0; j < len / 2; j++) {
      h_all_twiddles.push_back(w);
      w = ((int128)w * wlen_base) % reducer.mod;
    }
  }

  ll *d_all_twiddles;
  CUDA_CHECK(cudaMalloc((void **)&d_all_twiddles, total_twiddles * sizeof(ll)));
  CUDA_CHECK(cudaMemcpy(d_all_twiddles, h_all_twiddles.data(),
                        total_twiddles * sizeof(ll), cudaMemcpyHostToDevice));

  int num_blocks_twid =
      (total_twiddles + threads_per_block - 1) / threads_per_block;
  to_mont_kernel<Reducer><<<num_blocks_twid, threads_per_block>>>(
      d_all_twiddles, total_twiddles, reducer);
  CUDA_CHECK(cudaGetLastError());

  size_t twiddle_offset = 0;
  for (int len = 2; len <= n; len <<= 1) {
    ll *d_twiddles_stage = d_all_twiddles + twiddle_offset;

    int num_threads = n / 2;
    int num_blocks = (num_threads + threads_per_block - 1) / threads_per_block;

    ntt_kernel_optimized<Reducer><<<num_blocks, threads_per_block>>>(
        d_a, d_twiddles_stage, len, n, reducer);
    CUDA_CHECK(cudaGetLastError());

    twiddle_offset += len / 2;
  }

  CUDA_CHECK(cudaFree(d_all_twiddles));
}

std::vector<ll> multiply_ntt_gpu(std::vector<ll> &poly1, std::vector<ll> &poly2,
                                 ll mod, ll primitive_root,
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

  // GPU 版改为内核内部自行处理, 因此不再在 host 侧做 bit-reverse

  ll *d_p1, *d_p2;
  CUDA_CHECK(cudaMalloc((void **)&d_p1, n * sizeof(ll)));
  CUDA_CHECK(cudaMalloc((void **)&d_p2, n * sizeof(ll)));
  CUDA_CHECK(
      cudaMemcpy(d_p1, poly1.data(), n * sizeof(ll), cudaMemcpyHostToDevice));
  CUDA_CHECK(
      cudaMemcpy(d_p2, poly2.data(), n * sizeof(ll), cudaMemcpyHostToDevice));

  int threads = 256;
  int blocks = (n + threads - 1) / threads;

  if (method == "barrett") {
    BarrettReducer br(mod);
    ntt_gpu_optimized<BarrettReducer>(d_p1, n, false, br);
    ntt_gpu_optimized<BarrettReducer>(d_p2, n, false, br);
    pointwise_mult_kernel_optimized<BarrettReducer>
        <<<blocks, threads>>>(d_p1, d_p1, d_p2, n, br);
    ntt_gpu_optimized<BarrettReducer>(d_p1, n, true, br);
    ll n_inv = modInverse(n, mod);
    final_scaling_kernel_optimized<BarrettReducer>
        <<<blocks, threads>>>(d_p1, n, n_inv, br);
  } else if (method == "montgomery") {
    MontgomeryReducer mr(mod);

    to_mont_kernel<MontgomeryReducer><<<blocks, threads>>>(d_p1, n, mr);
    to_mont_kernel<MontgomeryReducer><<<blocks, threads>>>(d_p2, n, mr);

    ntt_gpu_montgomery(d_p1, n, false, mr);
    ntt_gpu_montgomery(d_p2, n, false, mr);

    pointwise_mult_kernel_optimized<MontgomeryReducer>
        <<<blocks, threads>>>(d_p1, d_p1, d_p2, n, mr);

    ntt_gpu_montgomery(d_p1, n, true, mr);

    ll n_inv = modInverse(n, mod);
    ll n_inv_mont = mr.to_mont(n_inv);
    final_scaling_kernel_optimized<MontgomeryReducer>
        <<<blocks, threads>>>(d_p1, n, n_inv_mont, mr);

    from_mont_kernel<MontgomeryReducer><<<blocks, threads>>>(d_p1, n, mr);

  } else {
    ntt_gpu_naive(d_p1, n, false, mod, primitive_root);
    ntt_gpu_naive(d_p2, n, false, mod, primitive_root);
    pointwise_mult_kernel<<<blocks, threads>>>(d_p1, d_p1, d_p2, n, mod);
    ntt_gpu_naive(d_p1, n, true, mod, primitive_root);
    ll n_inv = modInverse(n, mod);
    final_scaling_kernel<<<blocks, threads>>>(d_p1, n, n_inv, mod);
  }

  std::vector<ll> result(n);
  CUDA_CHECK(
      cudaMemcpy(result.data(), d_p1, n * sizeof(ll), cudaMemcpyDeviceToHost));

  CUDA_CHECK(cudaFree(d_p1));
  CUDA_CHECK(cudaFree(d_p2));

  result.resize(target_len);
  return result;
}