#include "values.h"

float compute_alpha(float cutoff, float sampling_freq) {
    return 1.0f/(float)std::tan((cutoff * std::numbers::pi) / sampling_freq);
}