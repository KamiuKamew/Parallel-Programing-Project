#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm> // For std::fill_n and std::min
#include <cassert>

// Assuming general/utils.h provides u64, u128, and expand_n
#include "../include/general/utils.h"
// The header with the function we want to test
#include "../include/pthread/ntt.h"

// Naive polynomial multiplication: res_full[k] = sum(a[i]*b[j]) where i+j=k
// Output coefficients are full precision (u128) before final modulo.
// res_full must be pre-allocated to at least n_coeffs_a + n_coeffs_b - 1,
// and zero-initialized if res_buffer_len is larger than that.
void poly_multiply_naive(const u64 *a, u64 n_coeffs_a,
                         const u64 *b, u64 n_coeffs_b,
                         u128 *res_full, u64 res_buffer_len)
{
    // Ensure res_full is zero-initialized up to res_buffer_len
    for (u64 k = 0; k < res_buffer_len; ++k)
    {
        res_full[k] = 0;
    }

    for (u64 i = 0; i < n_coeffs_a; ++i)
    {
        if (a[i] == 0)
            continue;
        for (u64 j = 0; j < n_coeffs_b; ++j)
        {
            if (b[j] == 0)
                continue;
            if (i + j < res_buffer_len)
            {
                res_full[i + j] += static_cast<u128>(a[i]) * b[j];
            }
        }
    }
}

// Helper to print polynomials
void print_poly(const char *name, const u64 *poly, u64 num_coeffs_to_print, u64 n_total_coeffs)
{
    std::cout << name << " (coeffs: ";
    for (u64 i = 0; i < std::min(num_coeffs_to_print, n_total_coeffs); ++i)
    {
        std::cout << poly[i] << " ";
    }
    if (num_coeffs_to_print < n_total_coeffs && n_total_coeffs > 0)
        std::cout << "...";
    if (n_total_coeffs == 0 && num_coeffs_to_print > 0)
        std::cout << "(empty)";
    else if (n_total_coeffs == 0 && num_coeffs_to_print == 0)
        std::cout << "(empty)";
    std::cout << ")" << std::endl;
}

int main()
{
    std::cout << "Testing pthread parallelized CRT NTT (poly_multiply_ntt_crt)..." << std::endl;
    std::cout << "===============================================================" << std::endl;

    // Test case 1: Small polynomials
    u64 a1_coeffs[] = {1, 2, 3};
    u64 b1_coeffs[] = {4, 5, 6};
    u64 n1_coeffs = sizeof(a1_coeffs) / sizeof(u64);
    assert(sizeof(b1_coeffs) / sizeof(u64) == n1_coeffs); // Assuming same length for simplicity in test
    u64 p1_final = CRT_MODS[0];                           // Use one of the CRT_MODS as final p: 998244353

    u64 result1_len_actual = 2 * n1_coeffs - 1;
    u64 n1_expanded = expand_n(result1_len_actual);

    u64 *res1_crt_pthread = new u64[n1_expanded];
    std::fill_n(res1_crt_pthread, n1_expanded, 0);

    u128 *res1_naive_full = new u128[n1_expanded];

    u64 *res1_naive_final = new u64[n1_expanded];
    std::fill_n(res1_naive_final, n1_expanded, 0);

    std::cout << std::endl
              << "Test Case 1: Small polynomials" << std::endl;
    std::cout << "Poly n_coeffs = " << n1_coeffs << ", Final Modulus p = " << p1_final << std::endl;
    std::cout << "Expected product coefficients = " << result1_len_actual << ", NTT expanded_len = " << n1_expanded << std::endl;
    print_poly("Poly A1", a1_coeffs, n1_coeffs, n1_coeffs);
    print_poly("Poly B1", b1_coeffs, n1_coeffs, n1_coeffs);

    poly_multiply_ntt_crt(a1_coeffs, b1_coeffs, res1_crt_pthread, n1_coeffs, p1_final);
    poly_multiply_naive(a1_coeffs, n1_coeffs, b1_coeffs, n1_coeffs, res1_naive_full, n1_expanded);
    for (u64 i = 0; i < n1_expanded; ++i)
    {
        res1_naive_final[i] = res1_naive_full[i] % p1_final;
    }

    print_poly("CRT_pthread Result1", res1_crt_pthread, std::min(result1_len_actual + 1, n1_expanded), n1_expanded);
    print_poly("Naive Result1", res1_naive_final, std::min(result1_len_actual + 1, n1_expanded), n1_expanded);

    bool test1_correct = true;
    for (u64 i = 0; i < n1_expanded; ++i)
    {
        if (res1_crt_pthread[i] != res1_naive_final[i])
        {
            std::cout << "ERROR (Test 1): Mismatch at index " << i
                      << ": CRT_Pthread = " << res1_crt_pthread[i]
                      << ", Naive = " << res1_naive_final[i] << std::endl;
            test1_correct = false;
        }
    }

    if (test1_correct)
    {
        std::cout << "Test Case 1 PASSED!" << std::endl;
    }
    else
    {
        std::cout << "Test Case 1 FAILED!" << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;

    delete[] res1_crt_pthread;
    delete[] res1_naive_full;
    delete[] res1_naive_final;

    // Test case 2: Slightly larger, different numbers
    u64 a2_coeffs[] = {123, 456, 789, 101, 112};
    u64 b2_coeffs[] = {234, 567, 890, 111, 222};
    u64 n2_coeffs = sizeof(a2_coeffs) / sizeof(u64);
    assert(sizeof(b2_coeffs) / sizeof(u64) == n2_coeffs);
    u64 p2_final = CRT_MODS[1]; // 1004535809

    u64 result2_len_actual = 2 * n2_coeffs - 1;
    u64 n2_expanded = expand_n(result2_len_actual);

    u64 *res2_crt_pthread = new u64[n2_expanded];
    std::fill_n(res2_crt_pthread, n2_expanded, 0);
    u128 *res2_naive_full = new u128[n2_expanded];
    u64 *res2_naive_final = new u64[n2_expanded];
    std::fill_n(res2_naive_final, n2_expanded, 0);

    std::cout << std::endl
              << "Test Case 2: Medium polynomials" << std::endl;
    std::cout << "Poly n_coeffs = " << n2_coeffs << ", Final Modulus p = " << p2_final << std::endl;
    std::cout << "Expected product coefficients = " << result2_len_actual << ", NTT expanded_len = " << n2_expanded << std::endl;
    print_poly("Poly A2", a2_coeffs, n2_coeffs, n2_coeffs);
    print_poly("Poly B2", b2_coeffs, n2_coeffs, n2_coeffs);

    poly_multiply_ntt_crt(a2_coeffs, b2_coeffs, res2_crt_pthread, n2_coeffs, p2_final);
    poly_multiply_naive(a2_coeffs, n2_coeffs, b2_coeffs, n2_coeffs, res2_naive_full, n2_expanded);
    for (u64 i = 0; i < n2_expanded; ++i)
    {
        res2_naive_final[i] = res2_naive_full[i] % p2_final;
    }

    print_poly("CRT_pthread Result2", res2_crt_pthread, std::min(result2_len_actual + 1, n2_expanded), n2_expanded);
    print_poly("Naive Result2", res2_naive_final, std::min(result2_len_actual + 1, n2_expanded), n2_expanded);

    bool test2_correct = true;
    for (u64 i = 0; i < n2_expanded; ++i)
    {
        if (res2_crt_pthread[i] != res2_naive_final[i])
        {
            std::cout << "ERROR (Test 2): Mismatch at index " << i
                      << ": CRT_Pthread = " << res2_crt_pthread[i]
                      << ", Naive = " << res2_naive_final[i] << std::endl;
            test2_correct = false;
        }
    }
    if (test2_correct)
    {
        std::cout << "Test Case 2 PASSED!" << std::endl;
    }
    else
    {
        std::cout << "Test Case 2 FAILED!" << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;

    delete[] res2_crt_pthread;
    delete[] res2_naive_full;
    delete[] res2_naive_final;

    // Test case 3: Edge case, n=1
    u64 a3_coeffs[] = {10};
    u64 b3_coeffs[] = {20};
    u64 n3_coeffs = sizeof(a3_coeffs) / sizeof(u64);
    assert(sizeof(b3_coeffs) / sizeof(u64) == n3_coeffs);
    u64 p3_final = 7340033; // A smaller prime, not necessarily in CRT_MODS from ntt.h

    u64 result3_len_actual = 2 * n3_coeffs - 1;
    u64 n3_expanded = expand_n(result3_len_actual);

    u64 *res3_crt_pthread = new u64[n3_expanded];
    std::fill_n(res3_crt_pthread, n3_expanded, 0);
    u128 *res3_naive_full = new u128[n3_expanded];
    u64 *res3_naive_final = new u64[n3_expanded];
    std::fill_n(res3_naive_final, n3_expanded, 0);

    std::cout << std::endl
              << "Test Case 3: n_coeffs = 1 (Edge case)" << std::endl;
    std::cout << "Poly n_coeffs = " << n3_coeffs << ", Final Modulus p = " << p3_final << std::endl;
    std::cout << "Expected product coefficients = " << result3_len_actual << ", NTT expanded_len = " << n3_expanded << std::endl;
    print_poly("Poly A3", a3_coeffs, n3_coeffs, n3_coeffs);
    print_poly("Poly B3", b3_coeffs, n3_coeffs, n3_coeffs);

    poly_multiply_ntt_crt(a3_coeffs, b3_coeffs, res3_crt_pthread, n3_coeffs, p3_final);
    poly_multiply_naive(a3_coeffs, n3_coeffs, b3_coeffs, n3_coeffs, res3_naive_full, n3_expanded);
    for (u64 i = 0; i < n3_expanded; ++i)
    {
        res3_naive_final[i] = res3_naive_full[i] % p3_final;
    }

    print_poly("CRT_pthread Result3", res3_crt_pthread, std::min(result3_len_actual + 1, n3_expanded), n3_expanded);
    print_poly("Naive Result3", res3_naive_final, std::min(result3_len_actual + 1, n3_expanded), n3_expanded);

    bool test3_correct = true;
    for (u64 i = 0; i < n3_expanded; ++i)
    {
        if (res3_crt_pthread[i] != res3_naive_final[i])
        {
            std::cout << "ERROR (Test 3): Mismatch at index " << i
                      << ": CRT_Pthread = " << res3_crt_pthread[i]
                      << ", Naive = " << res3_naive_final[i] << std::endl;
            test3_correct = false;
        }
    }
    if (test3_correct)
    {
        std::cout << "Test Case 3 PASSED!" << std::endl;
    }
    else
    {
        std::cout << "Test Case 3 FAILED!" << std::endl;
    }
    std::cout << "------------------------------------" << std::endl;

    delete[] res3_crt_pthread;
    delete[] res3_naive_full;
    delete[] res3_naive_final;

    std::cout << std::endl
              << "All tests for pthread CRT NTT completed." << std::endl;
    std::cout << "===============================================================" << std::endl;

    return 0;
}