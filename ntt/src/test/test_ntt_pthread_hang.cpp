#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
#include "../include/pthread_simple/ntt.h" // 调整路径
#include "../include/general/utils.h"      // 调整路径
#include "../include/general/op.h"         // 调整路径

// Helper function to print arrays for debugging
template <typename T>
void print_array(const std::string &name, T *arr, int n)
{
    std::cout << name << ": [";
    for (int i = 0; i < n; ++i)
    {
        std::cout << arr[i] << (i == n - 1 ? "" : ", ");
    }
    std::cout << "]" << std::endl;
}

// Helper function to compare arrays
template <typename T>
bool compare_arrays(T *arr1, T *arr2, int n)
{
    for (int i = 0; i < n; ++i)
    {
        if (arr1[i] != arr2[i])
        {
            std::cerr << "Mismatch at index " << i << ": " << arr1[i] << " != " << arr2[i] << std::endl;
            return false;
        }
    }
    return true;
}

int main()
{
    using T = int64_t; // Or your desired type

    T n_orig = 1 << 10;         // Original polynomial degree (e.g., 2^10)
    T p_mod = 998244353;        // A common NTT prime
    T omega_primitive_root = 3; // A primitive root for p_mod

    std::cout << "NTT Pthread Hang Test" << std::endl;
    std::cout << "Polynomial degree (original): " << n_orig << std::endl;
    std::cout << "Modulus P: " << p_mod << std::endl;
    std::cout << "Primitive root OMEGA: " << omega_primitive_root << std::endl;
    // std::cout << "Threads in pool: " << g_ntt_thread_pool.worker_pthread_ids.size() << std::endl; // Commented out as g_ntt_thread_pool is no longer global

    // For poly_multiply_ntt_pthread_simple
    T n_poly_mul = n_orig / 2; // Ensure 2*n_poly_mul-1 is reasonable
    std::vector<T> poly_a_orig(n_poly_mul);
    std::vector<T> poly_b_orig(n_poly_mul);
    // The output array ab must be large enough for n_expanded elements
    // n_expanded is expand_n(2 * n_poly_mul - 1)
    T temp_n_for_expansion = 2 * n_poly_mul - 1;
    if (temp_n_for_expansion < 1)
        temp_n_for_expansion = 1; // handle n_poly_mul = 0 or 1 edge cases for expansion
    T result_size_needed = expand_n(temp_n_for_expansion);

    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<T> distrib(0, p_mod - 1);

    for (T i = 0; i < n_poly_mul; ++i)
    {
        poly_a_orig[i] = distrib(rng);
        poly_b_orig[i] = distrib(rng);
    }

    const int NUM_ITERATIONS_POLY_MUL = 5;
    std::cout << "\nTesting poly_multiply_ntt_pthread_simple for " << NUM_ITERATIONS_POLY_MUL << " iterations..." << std::endl;
    for (int iter = 0; iter < NUM_ITERATIONS_POLY_MUL; ++iter)
    {
        std::cout << "Iteration " << iter + 1 << "/" << NUM_ITERATIONS_POLY_MUL << std::endl;

        // Prepare fresh copies for poly_multiply
        std::vector<T> a_current = poly_a_orig;
        std::vector<T> b_current = poly_b_orig;
        std::vector<T> ab_current(result_size_needed, 0); // Allocate n_expanded size

        poly_multiply_ntt_pthread_simple(a_current.data(), b_current.data(), ab_current.data(), n_poly_mul, p_mod, omega_primitive_root);

        // Basic check: just that it completes. A full correctness check would involve comparing to a naive poly multiply.
        std::cout << "poly_multiply_ntt_pthread_simple iteration " << iter + 1 << " completed." << std::endl;
    }
    std::cout << "poly_multiply_ntt_pthread_simple tests finished." << std::endl;

    // For direct ntt_forward_mont_pthread and ntt_inverse_mont_pthread
    // THIS SECTION IS NOW REMOVED AS THE FUNCTIONS ARE NO LONGER EXPOSED
    /*
    T n_expanded_direct = expand_n(n_orig);
    std::vector<T> data_orig(n_expanded_direct);
    std::vector<T> data_transformed(n_expanded_direct);
    std::vector<T> data_restored(n_expanded_direct);

    MontMod<T> montMod_direct(p_mod);
    T omega_mont_direct = montMod_direct.from_T(omega_primitive_root);
    T inv_omega_mont_direct = montMod_direct.inv(omega_mont_direct);

    for (T i = 0; i < n_expanded_direct; ++i)
    {
        data_orig[i] = distrib(rng);
    }

    const int NUM_ITERATIONS_DIRECT = 10;
    std::cout << "\\nTesting ntt_forward/inverse_mont_pthread for " << NUM_ITERATIONS_DIRECT << " iterations... (SKIPPED)" << std::endl;

    // The following loop is commented out because ntt_forward_mont_pthread and ntt_inverse_mont_pthread are no longer public.
    // for (int iter = 0; iter < NUM_ITERATIONS_DIRECT; ++iter)
    // {
    //     std::cout << "Iteration " << iter + 1 << "/" << NUM_ITERATIONS_DIRECT << std::endl;

    //     // Prepare data for transformation (in Montgomery form)
    //     for (T i = 0; i < n_expanded_direct; ++i)
    //     {
    //         data_transformed[i] = montMod_direct.from_T(data_orig[i]);
    //     }

    //     // Apply bit-reversal permutation before forward NTT
    //     bit_reverse_permute(data_transformed.data(), n_expanded_direct);

    //     // Forward NTT
    //     // ntt_forward_mont_pthread(data_transformed.data(), n_expanded_direct, p_mod, omega_mont_direct, montMod_direct);
    //     std::cout << "Forward NTT completed." << std::endl;

    //     // Inverse NTT
    //     // ntt_inverse_mont_pthread(data_transformed.data(), n_expanded_direct, p_mod, inv_omega_mont_direct, montMod_direct);
    //     std::cout << "Inverse NTT completed." << std::endl;

    //     // Data is now restored but still in Montgomery form and bit-reversed from inverse's perspective
    //     // Convert back from Montgomery form
    //     for (T i = 0; i < n_expanded_direct; ++i)
    //     {
    //         data_restored[i] = montMod_direct.to_T(data_transformed[i]);
    //     }

    //     // Apply bit-reversal permutation after inverse NTT
    //     bit_reverse_permute(data_restored.data(), n_expanded_direct);


    //     // Verification
    //     if (!compare_arrays(data_orig.data(), data_restored.data(), n_expanded_direct))
    //     {
    //         std::cerr << "Error: Data mismatch after forward and inverse NTT in iteration " << iter + 1 << std::endl;
    //         // print_array("Original", data_orig.data(), n_expanded_direct);
    //         // print_array("Restored", data_restored.data(), n_expanded_direct);
    //         return 1; // Indicate failure
    //     }
    //     std::cout << "Iteration " << iter + 1 << " verified successfully." << std::endl;
    // }
    */

    std::cout << "\nAll NTT Pthread hang tests completed successfully." << std::endl;
    return 0; // Indicate success
}