#include "../include/ntt.h"

#include <iostream>

int a[] = {1, 5, 5, 4, 0, 0, 0, 0};
int b[] = {4, 1, 5, 2, 0, 0, 0, 0};
int ab_naive[8], ab_ntt[8];

int main()
{
  poly_multiply_naive(a, b, ab_naive, 4, 998244353);
  poly_multiply_ntt(a, b, ab_ntt, 4, 998244353);

  std::cout << "naive: ";
  for (int i = 0; i < 8; ++i)
  {
    std::cout << ab_naive[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "ntt:   ";
  for (int i = 0; i < 8; ++i)
  {
    std::cout << ab_ntt[i] << " ";
  }
  std::cout << std::endl;

  for (int i = 0; i < 8; ++i)
  {
    if (ab_naive[i] != ab_ntt[i])
    {
      std::cout << "Test failed" << std::endl;
      return 1;
    }
  }
  std::cout << "Test passed" << std::endl;
}
