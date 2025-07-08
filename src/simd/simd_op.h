// This file is part of the Open Audio Live System project, a live audio environment
// Copyright (c) 2025 - Mathis DELGADO
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, version 3 of the License.
//
// This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

#ifndef OPENDSP_SIMD_OP_H
#define OPENDSP_SIMD_OP_H

#if defined(__x86_64__)
#include "immintrin.h"
#elif defined(__arm__) || defined(__aarch64__)
#include <arm_neon.h>
#endif

#include <ranges>
#include <cassert>

template<int __op_size>
float mulacc_no_simd(const std::ranges::subrange<float*>& a, const std::ranges::subrange<float*>& b) {
    float acc = 0.0f;

#pragma unroll
    for (int i = 0; i < __op_size; i++) {
        acc += a[i] * b[i];
    }

    return acc;
}

#if (defined(__arm__) || defined(__aarch64__))

// Finish mul acc and return
#define __INTERNAL_MULACC_RETF(op_a, op_b, zero) float32x4_t result = vfmaq_f32(zeros, op_a, op_b); return vaddvq_f32(result);

template<int __op_size>
float mulacc(const std::ranges::subrange<float*>& a, const std::ranges::subrange<float*>& b) {
    if constexpr (__op_size > 4) {
        return mulacc_no_simd<__op_size>(a, b);
    }

    return 0.0f;
}

template<>
inline float mulacc<1>(const std::ranges::subrange<float*>& a, const std::ranges::subrange<float*>& b) {
    // Trivial, no need for NEON
    return a[0] * b[0];
}

template<>
inline float mulacc<2>(const std::ranges::subrange<float*>& a, const std::ranges::subrange<float*>& b) {
    constexpr float32x4_t zeros = {0.0f, 0.0f, 0.0f, 0.0f};
    float32x4_t op_a = { a[0], a[1], 0.0f, 0.0f };
    float32x4_t op_b = { b[0], b[1], 0.0f, 0.0f };

    __INTERNAL_MULACC_RETF(op_a, op_b, zero);
}

template<>
inline float mulacc<3>(const std::ranges::subrange<float*>& a, const std::ranges::subrange<float*>& b) {
    constexpr float32x4_t zeros = {0.0f, 0.0f, 0.0f, 0.0f};
    float32x4_t op_a = { a[0], a[1], a[2], 0.0f };
    float32x4_t op_b = { b[0], b[1], b[2], 0.0f };

    __INTERNAL_MULACC_RETF(op_a, op_b, zero);
}

template<>
inline float mulacc<4>(const std::ranges::subrange<float*>& a, const std::ranges::subrange<float*>& b) {
    constexpr float32x4_t zeros = {0.0f, 0.0f, 0.0f, 0.0f};
    float32x4_t op_a = { a[0], a[1], a[2], a[3] };
    float32x4_t op_b = { b[0], b[1], b[2], b[3] };

    __INTERNAL_MULACC_RETF(op_a, op_b, zero);
}

#else

template<int __op_size>
float mulacc(const std::ranges::subrange<float*>& a, const std::ranges::subrange<float*>& b) {
    return mulacc_no_simd(a, b);
}

#endif

#endif //OPENDSP_SIMD_OP_H
