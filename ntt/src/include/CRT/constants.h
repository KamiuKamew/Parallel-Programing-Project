#pragma once

#include "../general/type.h"

// CRT constants
// Using static to keep them file-local. Names changed to avoid potential clashes if included elsewhere.
static const u32 CRT_MODS[] = {998244353, 1004535809}; // P0, P1
static const u32 CRT_ROOTS[] = {3, 3};                 // G0, G1
static const int CRT_NUMS = sizeof(CRT_MODS) / sizeof(CRT_MODS[0]);
