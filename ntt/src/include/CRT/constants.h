#pragma once

#include "../general/type.h"

const int ROOT = 3;
// NTT素数，这些素数的原根均为 ROOT (3)
// 来自 https://blog.miskcoo.com/2014/07/fft-prime-table
const u64 NTT_PRIMES[] = {167772161, 469762049, 998244353, 1004535809};
const int NUM_NTT_PRIMES = sizeof(NTT_PRIMES) / sizeof(NTT_PRIMES[0]);
