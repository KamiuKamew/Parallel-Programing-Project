#pragma once

int mod_add(int a, int b, int mod) { return (a + b) % mod; }
int mod_sub(int a, int b, int mod) { return (a - b + mod) % mod; }
int mod_mul(int a, int b, int mod) { return (1LL * a * b) % mod; }
int mod_pow(int base, int exp, int mod) {
  int result = 1;
  while (exp > 0) {
    if (exp & 1) {
      result = mod_mul(result, base, mod);
    }
    base = mod_mul(base, base, mod);
    exp >>= 1;
  }
  return result;
}
int mod_inv(int x, int mod) { return mod_pow(x, mod - 2, mod); }
