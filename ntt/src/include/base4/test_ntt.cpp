#include "ntt-all-in-one.h"
#include <iostream>
#include <vector>
#include <numeric>

// Helper function to print arrays
template <typename T>
void print_array(const char* name, T* arr, T n) {
    std::cout << name << " = [ ";
    for (T i = 0; i < n; ++i) {
        std::cout << arr[i] << (i == n - 1 ? "" : ", ");
    }
    std::cout << " ]\n";
}

// Test Mod class operations
template <typename T>
void test_mod_operations(T p) {
    Mod<T> mod(p);
    std::cout << "\nTesting Mod<" << typeid(T).name() << "> with p = " << p << "\n";
    T a = p / 2 + 1;
    T b = p / 3 + 1;
    std::cout << "a = " << a << ", b = " << b << "\n";
    std::cout << "add(a, b) = " << mod.add(a, b) << " (expected: " << (a + b) % p << ")\n";
    std::cout << "sub(a, b) = " << mod.sub(a, b) << " (expected: " << (a >= b ? a - b : a + p - b) << ")\n";
    std::cout << "mul(a, b) = " << mod.mul(a, b) << " (expected: " << ((u64)a * b) % p << ")\n";
    std::cout << "pow(a, 3) = " << mod.pow(a, 3) << " (expected: " << mod.mul(mod.mul(a, a), a) << ")\n";
    std::cout << "inv(a) = " << mod.inv(a) << " (expected: " << mod.pow(a, p - 2) << ")\n";
    std::cout << "reduce(a + p) = " << mod.reduce(a + p) << " (expected: " << a % p << ")\n";
}

// Test NTT forward and inverse for a given type T
template <typename T>
void test_ntt_transform(T p, T omega, const std::vector<T>& input_poly, const char* type_name, bool radix4 = false) {
    std::cout << "\nTesting " << (radix4 ? "Radix-4 " : "") << "NTT for " << type_name << " with p = " << p << ", omega = " << omega << "\n";
    T n = input_poly.size();
    T n_expanded = expand_n(n);
    std::vector<T> a_expanded_vec(n_expanded);
    for (T i = 0; i < n; ++i) a_expanded_vec[i] = input_poly[i];
    for (T i = n; i < n_expanded; ++i) a_expanded_vec[i] = 0;

    T* a_expanded = a_expanded_vec.data();

    std::vector<T> original_a(a_expanded_vec);

    std::cout << "Original polynomial: ";
    print_array("", original_a.data(), n_expanded);

    // Forward NTT
    bit_reverse_permute(a_expanded, n_expanded);
    if (radix4) {
        ntt_forward_radix4(a_expanded, n_expanded, p, omega);
    } else {
        ntt_forward(a_expanded, n_expanded, p, omega);
    }
    std::cout << "Forward NTT result: ";
    print_array("", a_expanded, n_expanded);

    // Inverse NTT
    if (radix4) {
        ntt_inverse_radix4(a_expanded, n_expanded, p, omega);
    } else {
        ntt_inverse(a_expanded, n_expanded, p, omega);
    }
    std::cout << "Inverse NTT result: ";
    print_array("", a_expanded, n_expanded);

    // Verify if inverse result matches original
    bool success = true;
    for (T i = 0; i < n_expanded; ++i) {
        if (a_expanded[i] != original_a[i]) {
            success = false;
            break;
        }
    }
    std::cout << "NTT transform and inverse: " << (success ? "SUCCESS" : "FAILED") << "\n";
}

// Test polynomial multiplication using NTT
template <typename T>
void test_poly_multiply(T p, T omega, const std::vector<T>& a_poly, const std::vector<T>& b_poly, const char* type_name) {
    std::cout << "\nTesting poly_multiply_ntt for " << type_name << " with p = " << p << ", omega = " << omega << "\n";

    T n_a = a_poly.size();
    T n_b = b_poly.size();
    T n_result = n_a + n_b - 1;

    std::vector<T> a_vec(a_poly);
    std::vector<T> b_vec(b_poly);
    std::vector<T> ab_vec(n_result);

    print_array("a", a_vec.data(), n_a);
    print_array("b", b_vec.data(), n_b);

    poly_multiply_ntt(a_vec.data(), b_vec.data(), ab_vec.data(), n_a, n_b, p, omega);

    std::cout << "NTT multiplication result: ";
    print_array("ab", ab_vec.data(), n_result);

    // Manual multiplication for verification
    Mod<T> mod(p);
    std::vector<T> expected_ab(n_result, 0);
    for (T i = 0; i < n_a; ++i) {
        for (T j = 0; j < n_b; ++j) {
            expected_ab[i + j] = mod.add(expected_ab[i + j], mod.mul(a_vec[i], b_vec[j]));
        }
    }
    std::cout << "Expected multiplication result: ";
    print_array("expected_ab", expected_ab.data(), n_result);

    bool success = true;
    for (T i = 0; i < n_result; ++i) {
        if (ab_vec[i] != expected_ab[i]) {
            success = false;
            break;
        }
    }
    std::cout << "Polynomial multiplication: " << (success ? "SUCCESS" : "FAILED") << "\n";
}

int main() {
    // Test u32
    test_mod_operations<u32>(998244353);
    test_ntt_transform<u32>(998244353, 3, {1, 2, 3, 4}, "u32");
    test_poly_multiply<u32>(998244353, 3, {1, 2}, {3, 4}, "u32"); // (1+2x)(3+4x) = 3 + 10x + 8x^2
    test_poly_multiply<u32>(998244353, 3, {1, 2, 3}, {4, 5, 6}, "u32"); // (1+2x+3x^2)(4+5x+6x^2) = 4 + 13x + 28x^2 + 27x^3 + 18x^4

    // Test Radix-4 NTT for u32
    // Modulus for Radix-4 NTT must satisfy p = 1 (mod 4)
    // 998244353 = 1 (mod 4) is true. So we can use it.
    test_ntt_transform<u32>(998244353, 3, {1, 2, 3, 4}, "u32", true);
    test_ntt_transform<u32>(998244353, 3, {1, 2, 3, 4, 5, 6, 7, 8}, "u32", true);

    return 0;
}


