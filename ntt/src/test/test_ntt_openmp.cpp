#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <chrono>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <cmath>  // For std::isspace in C++98 mode, or for general use
#include <cctype> // For std::isspace with a specific locale if needed

// Assuming basicMultiThread.h is in ntt/src/include/thread/
// Relative path from ntt/src/test/ to ntt/src/include/thread/ is ../include/thread/
#include "../include/thread/basicMultiThread.h"

// Function to read a polynomial from an input stream
bool read_poly(std::ifstream &ifs, int n, std::vector<int> &poly)
{
    poly.resize(n);
    for (int i = 0; i < n; ++i)
    {
        if (!(ifs >> poly[i]))
        {
            // std::cerr << "Error reading coefficient " << i << " for polynomial of size " << n << std::endl;
            return false;
        }
    }
    return true;
}

// Function to read expected output coefficients
bool read_expected_output(std::ifstream &ifs, int count, std::vector<int> &expected_poly)
{
    expected_poly.resize(count);
    for (int i = 0; i < count; ++i)
    {
        if (!(ifs >> expected_poly[i]))
        {
            // std::cerr << "Error reading expected coefficient " << i << " for count " << count << std::endl;
            return false;
        }
    }
    return true;
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file.in> <output_file.out>" << std::endl;
        return 1; // General error
    }

    std::string input_filename = argv[1];
    std::string output_filename = argv[2];

    std::ifstream ifs(input_filename);
    if (!ifs.is_open())
    {
        std::cerr << "Error: Could not open input file " << input_filename << std::endl;
        return 2; // Specific error code for input file open error
    }

    int n_input;
    u32 p_modulus; // Changed to u32 to match MontMod constructor and poly_multiply_ntt
    if (!(ifs >> n_input >> p_modulus))
    {
        std::cerr << "Error: Could not read n and p from " << input_filename << std::endl;
        ifs.close();
        return 3; // Specific error code for n, p read error
    }

    std::vector<int> poly_a_int, poly_b_int;
    if (!read_poly(ifs, n_input, poly_a_int) || !read_poly(ifs, n_input, poly_b_int))
    {
        std::cerr << "Error: Could not read polynomials from " << input_filename << std::endl;
        ifs.close();
        return 4; // Specific error code for polynomial read error
    }
    ifs.close();

    u32 actual_n_for_poly_multiply = n_input;
    if (n_input <= 0)
    {
        // std::cerr << "Warning: n_input is " << n_input << ". Adjusted to 0 for poly_multiply_ntt length parameter." << std::endl;
        actual_n_for_poly_multiply = 0;
    }

    u32 n_expanded;
    int expected_coeffs_count;

    if (actual_n_for_poly_multiply == 0)
    {
        n_expanded = 0;
        expected_coeffs_count = 0;
    }
    else
    {
        n_expanded = expand_n(2 * actual_n_for_poly_multiply - 1);
        expected_coeffs_count = 2 * actual_n_for_poly_multiply - 1;
    }
    if (n_expanded == 0 && actual_n_for_poly_multiply > 0)
    { // e.g. n_input=1, 2*1-1=1, expand_n(1)=1
        n_expanded = 1;
    }

    std::vector<int> result_coeffs_int(n_expanded > 0 ? n_expanded : 1, 0); // Ensure min size 1

    auto start_time = std::chrono::high_resolution_clock::now();
    if (actual_n_for_poly_multiply > 0)
    {
        poly_multiply_ntt(poly_a_int.data(), poly_b_int.data(), result_coeffs_int.data(), actual_n_for_poly_multiply, p_modulus);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration_us = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

    std::ifstream ofs_expected(output_filename);
    if (!ofs_expected.is_open())
    {
        std::cerr << "Error: Could not open expected output file " << output_filename << std::endl;
        return 5; // Specific error code for output file open error
    }

    std::vector<int> expected_poly;
    bool read_out_ok = true;
    if (expected_coeffs_count > 0)
    {
        if (!read_expected_output(ofs_expected, expected_coeffs_count, expected_poly))
        {
            std::cerr << "Error: Could not read expected output from " << output_filename << std::endl;
            read_out_ok = false;
        }
    }
    else
    {
        char ch_test;
        bool has_non_whitespace = false;
        while (ofs_expected.get(ch_test))
        {
            if (!std::isspace(static_cast<unsigned char>(ch_test)))
            {
                has_non_whitespace = true;
                break;
            }
        }
        if (has_non_whitespace)
        {
            std::cerr << "Warning: Expected 0 coefficients, but output file " << output_filename << " has non-whitespace content." << std::endl;
        }
    }
    ofs_expected.close();
    if (!read_out_ok && expected_coeffs_count > 0)
    {
        return 6; // Specific error code for expected output read error
    }

    bool pass = true;
    if (expected_coeffs_count > 0)
    {
        if (result_coeffs_int.size() < (size_t)expected_coeffs_count)
        {
            pass = false;
            std::cerr << "Error: Result vector too small. Got " << result_coeffs_int.size() << ", expected " << expected_coeffs_count << std::endl;
        }
        else
        {
            for (int i = 0; i < expected_coeffs_count; ++i)
            {
                if (result_coeffs_int[i] != expected_poly[i])
                {
                    pass = false;
                    break;
                }
            }
        }
    }
    else
    {
        for (size_t i = 0; i < result_coeffs_int.size(); ++i)
        {
            if (result_coeffs_int[i] != 0)
            { // If expecting 0 coeffs, any non-zero in result (even if padded) is a fail.
                pass = false;
                break;
            }
        }
    }

    if (!pass)
    {
        std::cerr << "FAIL: " << input_filename << std::endl;
        if (expected_coeffs_count > 0)
        {
            std::cerr << "  Expected (" << expected_coeffs_count << "): ";
            for (int i = 0; i < expected_coeffs_count; ++i)
                std::cerr << expected_poly[i] << (i == expected_coeffs_count - 1 ? "" : " ");
            std::cerr << std::endl;
            std::cerr << "  Got      (" << std::min((size_t)expected_coeffs_count, result_coeffs_int.size()) << "): ";
            for (size_t i = 0; i < std::min((size_t)expected_coeffs_count, result_coeffs_int.size()); ++i)
                std::cerr << result_coeffs_int[i] << (i == std::min((size_t)expected_coeffs_count, result_coeffs_int.size()) - 1 ? "" : " ");
            std::cerr << std::endl;
            if (result_coeffs_int.size() > (size_t)expected_coeffs_count && (size_t)expected_coeffs_count > 0)
            {
                std::cerr << "  (Result had more coeffs, up to " << result_coeffs_int.size() << ")" << std::endl;
            }
        }
        else
        {
            std::cerr << "  Expected: (0 coefficients or all zeros)" << std::endl;
            std::cerr << "  Got (size " << result_coeffs_int.size() << "): ";
            for (size_t i = 0; i < std::min((size_t)5, result_coeffs_int.size()); ++i)
                std::cerr << result_coeffs_int[i] << " ";
            if (result_coeffs_int.size() > 5)
                std::cerr << "...";
            std::cerr << std::endl;
        }
        return 7; // Specific error code for correctness failure
    }

    // If PASS, print only the duration in microseconds to stdout
    std::cout << duration_us.count() << std::endl;
    return 0; // Success
}