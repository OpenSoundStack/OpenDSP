// This file is part of the Open Audio Live System project, a live audio environment
// Copyright (c) 2026 - Mathis DELGADO
//
// This project is distributed under the Creative Commons CC-BY-NC-SA licence. https://creativecommons.org/licenses/by-nc-sa/4.0

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
            float factor = 0.5f;

            if ((time > low_phase_dur_ms) && (time < (low_phase_dur_ms + high_phase_dur_ms))) {
                factor = 1.0f;
            }

            comp_test.push_back(
                sin(2.0f * 3.141592f * 5000.0f * i * (1.0f/96000.0f)) * factor
            );
        }
    }

    int attack = 4;
    Enveloppe env_in{10, 96000};
    std::vector<float> enveloppe_in;

    for (auto& s : comp_test) {
        enveloppe_in.push_back(10 * log10(env_in.push_sample(s)));
    }

    std::vector<float> gain_red;
    std::vector<float> signal_out;
    Dynamics dyn{[](float level_db) {
        float threshold = -4;
        float ratio = 2.0f;

        if (level_db > threshold) {
            return -(level_db - threshold) * ratio;
        } else {
            return 0.0f;
        }
    }, attack, 70, 50, 96000};

    for (auto& s : comp_test) {
        gain_red.push_back(20.0f * log10(dyn.push_sample(s)));
    }

    for (int i = 0; i < comp_test.size(); i++) {
        signal_out.push_back(comp_test[i] * std::pow(10, gain_red[i] / 20));
    }

    plt::plot(gain_red);
   // plt::plot(enveloppe_in);
    plt::show();

    plt::plot(signal_out);
    plt::show();

    return 0;
}