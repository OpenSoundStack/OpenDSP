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

#include <iostream>
#include <chrono>
#include <cmath>

#include "dynamics/enveloppe.h"
#include "dynamics/dynamics.h"

#define WITHOUT_NUMPY
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

using namespace std::chrono;
using namespace std::chrono_literals;

int main() {
    std::vector<float> comp_test{};

    // Generating a pulse for dynamics analysis
    int low_phase_dur_ms = 150;
    int high_phase_dur_ms = 60;
    int total_dur_samples = (2 * low_phase_dur_ms + high_phase_dur_ms) * 96;

    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < total_dur_samples; i++) {
            int time = i / 96;
            float factor = 1.0f;

            if ((time > low_phase_dur_ms) && (time < (low_phase_dur_ms + high_phase_dur_ms))) {
                factor = 2.0f;
            }

            comp_test.push_back(
                sin(2.0f * 3.141592f * 5000.0f * i * (1.0f/96000.0f)) * factor
            );
        }
    }

    int attack = 70;
    Enveloppe env_in{attack, 96000};
    std::vector<float> enveloppe_in;

    for (auto& s : comp_test) {
        enveloppe_in.push_back((env_in.push_sample(s)));
    }

    std::vector<float> gain_red;
    std::vector<float> signal_out;
    Dynamics dyn{[](float level_lin) {
        float threshold = 0.72f;
        if (level_lin > threshold) {
            return 0.5f * level_lin + 0.5f * threshold;
        } else {
            return 1.0f * level_lin;
        }
    }, attack, 150, 10, 96000};

    for (auto& s : comp_test) {
        gain_red.push_back((dyn.push_sample(s)));
    }

    for (int i = 0; i < comp_test.size(); i++) {
        signal_out.push_back(comp_test[i] * gain_red[i]);
    }

    plt::plot(gain_red);
    //plt::plot(enveloppe_out);
    plt::plot(enveloppe_in);
    plt::show();

    plt::plot(signal_out);
    plt::show();

    return 0;
}