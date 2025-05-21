#include "../include/ntt.h"
// #include "../include/simd/ntt.h"
#include "../include/CRT/ntt.h"

#include <iostream>

u64 a[] = {1, 5, 5, 4, 0, 0, 0, 0};
u64 b[] = {4, 1, 5, 5, 0, 0, 0, 0};
u64 ab_naive[8], ab_ntt[8], ab_ntt_simd[8], ab_ntt_crt[8];

int main()
{
  // poly_multiply_naive(a, b, ab_naive, 4, 998244353);
  poly_multiply_ntt<u64>(a, b, ab_ntt, 4, 998244353);
  // poly_multiply_ntt_simd(a, b, ab_ntt_simd, 4, 998244353);
  poly_multiply_ntt_crt(a, b, ab_ntt_crt, 4, 998244353);

  // std::cout << "naive:   ";
  // for (int i = 0; i < 8; ++i)
  // {
  //   std::cout << ab_naive[i] << " ";
  // }
  // std::cout << std::endl;

  std::cout << "ntt:     ";
  for (int i = 0; i < 8; ++i)
  {
    std::cout << ab_ntt[i] << " ";
  }
  std::cout << std::endl;

  // std::cout << "ntt_simd:";
  // for (int i = 0; i < 8; ++i)
  // {
  //   std::cout << ab_ntt_simd[i] << " ";
  // }
  // std::cout << std::endl;

  std::cout << "ntt_crt:";
  for (int i = 0; i < 8; ++i)
  {
    std::cout << ab_ntt_crt[i] << " ";
  }
  std::cout << std::endl;
}
