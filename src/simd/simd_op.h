#ifndef OPENDSP_SIMD_OP_H
#define OPENDSP_SIMD_OP_H

#if defined(__x86_64__)
#include "immintrin.h"
#elif defined(__arm__) || defined(__aarch64__)
#include <arm_neon.h>
#endif

#endif //OPENDSP_SIMD_OP_H
