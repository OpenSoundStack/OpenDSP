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

#ifndef OPENDSP_BANDPASS_H
#define OPENDSP_BANDPASS_H

#include "filter/iirfilter.h"
#include "utils/values.h"

class BPF_2ord {
public:
    BPF_2ord(float f0, float Q, float sampling_freq) {
        m_f0 = f0;
        m_Q = Q;
        m_sampling_freq = sampling_freq;

        m_alpha = compute_alpha(f0, sampling_freq);

        init_filter();
    }

    ~BPF_2ord() = default;

    float push_sample(const float& s) {
        return m_filter.push_sample(s);
    }

    float get_output() const {
        return m_filter.get_output();
    }

    void set_central_freq(float f0) {
        m_f0 = f0;
        m_alpha = compute_alpha(m_f0, m_sampling_freq);

        update_filter();
    }

    void set_quality_factor(float Q) {
        m_Q = Q;
        update_filter();
    }

    IIRFilter<2>& get_filter() {
        return m_filter;
    }

private:
    void init_filter() {
        auto weights = compute_weights();
        m_filter = IIRFilter<2>(weights[0], weights[1]);
    }

    void update_filter() {
        auto weights = compute_weights();
        m_filter.set_weights(weights[0], weights[1]);
    }

    std::array<std::array<float, 3>, 2> compute_weights() {
        float alpha2 = m_alpha * m_alpha;
        float aQ = m_alpha / m_Q;
        float inv_com_den = 1.0f + aQ + alpha2; // Inverse common denominator
        inv_com_den = 1.0f / inv_com_den;

        std::array<float, 3> xweights = {
                aQ * inv_com_den,
                0.0f,
                (-aQ) * inv_com_den
        };

        std::array<float, 3> yweights = {
                1.0f,
                (2.0f - (2.0f * alpha2)) * inv_com_den,
                (1.0f - aQ + alpha2) * inv_com_den
        };

        return {xweights, yweights};
    }

    float m_f0;
    float m_Q;
    float m_sampling_freq;
    float m_alpha;

    IIRFilter<2> m_filter;
};

#endif //OPENDSP_BANDPASS_H
