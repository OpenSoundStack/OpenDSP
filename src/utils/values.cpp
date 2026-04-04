// This file is part of the Open Audio Live System project, a live audio environment
// Copyright (c) 2026 - Mathis DELGADO
//
// This project is distributed under the Creative Commons CC-BY-NC-SA licence. https://creativecommons.org/licenses/by-nc-sa/4.0

#include "values.h"

float compute_alpha(float cutoff, float sampling_freq) {
    return 1.0f/(float)std::tan((cutoff * std::numbers::pi) / sampling_freq);
}