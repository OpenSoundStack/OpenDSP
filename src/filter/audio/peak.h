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

#ifndef OPENDSP_PEAK_H
#define OPENDSP_PEAK_H

#include "filter/iirfilter.h"
#include "utils/values.h"

class PeakFilter {
public:
    PeakFilter(float fc, float Q, float gain_db, float sampling_freq) {
        m_fc = fc;
        m_Q = Q;
        m_gain = gain_db;
        m_sampling_freq = sampling_freq;

        init_filter();
    }

    ~PeakFilter() = default;

    float push_sample(const float& s) {
        return m_filter.push_sample(s);
    }

    float get_output() const {
        return m_filter.get_output();
    }

    void set_cutoff(float fc) {
        m_fc = fc;
        update_filter();
    }

    void set_quality_factor(float Q) {
        m_Q = Q;
        update_filter();
    }

    void set_gain(float gain) {
        m_gain = gain;
        update_filter();
    }

    float get_gain() {
        return m_gain;
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
        float A = std::pow(10.0f, m_gain / 40.0f);
        float rQ = m_Q;
        float w0 = 2.0f * std::numbers::pi * (m_fc / m_sampling_freq);
        float alpha = std::sin(w0) / (2.0f * rQ);

        float inv_com_den = 1.0f + (alpha / A);
        inv_com_den = 1.0f / inv_com_den;

        float B = (-2.0f * std::cos(w0)) * inv_com_den;
        std::array<float, 3> xweights = {
                (1.0f + (alpha * A)) * inv_com_den,
                B,
                (1.0f - (alpha * A)) * inv_com_den
        };

        std::array<float, 3> yweights = {
                1.0f,
                B,
                (1.0f - (alpha / A)) * inv_com_den
        };

        return {xweights, yweights};
    }

    float m_Q; // Quality factor
    float m_fc; // Cutoff frequency
    float m_gain;
    float m_sampling_freq;

    IIRFilter<2> m_filter;
};

#endif //OPENDSP_PEAK_H
