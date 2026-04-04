// This file is part of the Open Audio Live System project, a live audio environment
// Copyright (c) 2026 - Mathis DELGADO
//
// This project is distributed under the Creative Commons CC-BY-NC-SA licence. https://creativecommons.org/licenses/by-nc-sa/4.0

#ifndef MATMUL_H
#define MATMUL_H

#include "simd_op.h"

struct Mat4x4 {
    float coefs[4][4];
};

struct Vec4 {
    float elems[4];
};

inline float vec4_product(const Vec4& v) {
    float acc = v.elems[0];
    for (int i = 1; i < 4; i++) {
        acc *= v.elems[i];
    }

    return acc;
}

#if defined(SIMD_ENABLE)

inline Vec4 mat4_vec4_mul(const Mat4x4& m, const Vec4& v) {
    Vec4 outvec = {};

    auto mat = vld1q_f32_x4((float*)m.coefs);
    float32x4_t vec = vld1q_f32(v.elems);
    float32x4_t zero = {0, 0, 0, 0};

    for (int i = 0; i < 4; i++) {
        float32x4_t result = vfmaq_f32(zero, vec, mat.val[i]);
        outvec.elems[i] = vaddvq_f32(result);
    }

    return outvec;
}

inline Vec4 mat4_vec4_mul(const Mat4x4& m, float* v) {
    Vec4 outvec = {};

    auto mat = vld1q_f32_x4((float*)m.coefs);
    float32x4_t vec = vld1q_f32(v);
    float32x4_t zero = {0, 0, 0, 0};

    for (int i = 0; i < 4; i++) {
        float32x4_t result = vfmaq_f32(zero, vec, mat.val[i]);
        outvec.elems[i] = vaddvq_f32(result);
    }

    return outvec;
}

inline Vec4 vec4_sum(const Vec4& v1, const Vec4& v2) {
    Vec4 outvec = {};

    float32x4_t neon_v1 = vld1q_f32(v1.elems);
    float32x4_t neon_v2 = vld1q_f32(v2.elems);

    float32x4_t res = vaddq_f32(neon_v1, neon_v2);
    vst1q_f32(outvec.elems, res);

    return outvec;
}

inline Vec4 vec4_diff(const Vec4& v1, const Vec4& v2) {
    Vec4 outvec = {};

    float32x4_t neon_v1 = vld1q_f32(v1.elems);
    float32x4_t neon_v2 = vld1q_f32(v2.elems);

    float32x4_t res = vsubq_f32(neon_v1, neon_v2);
    vst1q_f32(outvec.elems, res);

    return outvec;
}

#else

inline Vec4 vec4_sum(const Vec4& v1, const Vec4& v2) {
    Vec4 outvec = {};

    for (int i = 0; i < 4; i++) {
        outvec.elems[i] = v1.elems[i] + v2.elems[i];
    }

    return outvec;
}

inline Vec4 vec4_diff(const Vec4& v1, const Vec4& v2) {
    Vec4 outvec = {};

    for (int i = 0; i < 4; i++) {
        outvec.elems[i] = v1.elems[i] - v2.elems[i];
    }

    return outvec;
}

inline Vec4 mat4_vec4_mul(const Mat4x4& m, const Vec4& v) {
    Vec4 outvec = {};

    for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 4; col++) {
            outvec.elems[row] += m.coefs[row][col] * v.elems[col];
        }
    }

    return outvec;
}

inline Vec4 mat4_vec4_mul(const Mat4x4& m, float* v) {
    Vec4 outvec = {};

    for (int row = 0; row < 4; row++) {
        for (int col = 0; col < 4; col++) {
            outvec.elems[row] += m.coefs[row][col] * v[col];
        }
    }

    return outvec;
}

#endif

#endif //MATMUL_H
