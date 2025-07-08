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

template<int op_size__, int init__>
float mulacc_no_simd(const float* a, const float* b) {
    float acc = 0.0f;
    for (int i = init__; i < op_size__; i++) {
        acc += a[i] * b[i];
    }

    return acc;
}

#if (defined(__arm__) || defined(__aarch64__))

// Finish mul acc and return
#define __INTERNAL_MULACC_RETF(op_a, op_b, zero) float32x4_t result = vfmaq_f32(zeros, op_a, op_b); return vaddvq_f32(result);

template<int __op_size, int init__>
float mulacc(const float* a, const float* b) {
    if (__op_size > 4) {
        assert(false && "Not implemented yet");
        return 0.0f;
    } else {
        float32x4_t zeros = {0.0f, 0.0f, 0.0f, 0.0f};

        if constexpr (__op_size - init__ == 1) {
            // Trivial, no need for NEON
            return a[init__] * b[init__];
        } else if (__op_size - init__ == 2) {
            float32x4_t op_a = { a[init__], a[init__ + 1], 0.0f, 0.0f };
            float32x4_t op_b = { b[init__], b[init__ + 1], 0.0f, 0.0f };

            __INTERNAL_MULACC_RETF(op_a, op_b, zero);
        } else if (__op_size - init__ == 3) {
            float32x4_t op_a = { a[init__], a[init__ + 1], a[init__ + 2], 0.0f};
            float32x4_t op_b = { b[init__], b[init__ + 1], b[init__ + 2], 0.0f};

            __INTERNAL_MULACC_RETF(op_a, op_b, zero);
        } else if (__op_size - init__ == 4) {
            float32x4_t op_a = { a[init__], a[init__ + 1], a[init__ + 2], a[init__ + 3]};
            float32x4_t op_b = { b[init__], b[init__ + 1], b[init__ + 2], b[init__ + 3]};

            __INTERNAL_MULACC_RETF(op_a, op_b, zero);
        }

    }

    return 0.0f;
}

#else

template<int __op_size>
float mulacc(const std::ranges::subrange<float*>& a, const std::ranges::subrange<float*>& b) {
    return mulacc_no_simd<__op_size>(a, b);
}

#endif

#endif //OPENDSP_SIMD_OP_H
