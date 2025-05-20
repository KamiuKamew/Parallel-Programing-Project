#pragma once

#include <stdint.h>
#include <type_traits>

using u32 = uint32_t;
using u64 = uint64_t;
using u128 = __uint128_t;

using u32_mont = u32;
using u64_mont = u64;

using t_widen_default = u128;
template <typename T>
using t_widen = typename std::conditional<std::is_same<T, u32>::value, u64,
                                          typename std::conditional<std::is_same<T, u64>::value, u128, t_widen_default>::type>::type;
