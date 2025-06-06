#pragma once

#include "../general/op.h"

template <typename T>
class BarrettMod : public Mod<T>
{
    using T2 = t_widen<T>; // TODO：当 T = u128 时，a * b 可能溢出

public:
    BarrettMod(T _mod, int iterations = 0) : Mod<T>(_mod)
    {
        barrett_k = 64;

        T2 r = (T2(1) << barrett_k) / _mod;
        for (int i = 0; i < iterations; ++i)
        {
            T2 numerator = (T2(1) << barrett_k) - r * _mod;
            T2 correction = numerator / _mod;
            r += correction;
        }

        barrett_r = r;
    }

    T mul(T a, T b) const
    {
        T2 z = T2(a) * b;
        T2 q = (z * barrett_r) >> barrett_k;
        T res = T(z - q * this->mod);
        if (res >= this->mod)
            res -= this->mod;
        if (res >= this->mod)
            res -= this->mod; // 再减一次以确保结果 < mod
        return res;
    }

private:
    T2 barrett_r;
    int barrett_k;
};
