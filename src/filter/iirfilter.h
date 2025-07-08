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

#ifndef OPENDSP_IIRFILTER_H
#define OPENDSP_IIRFILTER_H

#include <array>
#include <numbers>

#include "utils/sample_buffer.h"
#include "utils/values.h"
#include "simd/simd_op.h"

template<int order__>
class IIRFilter {
public:
    IIRFilter(const std::array<float, order__ + 1>& xweights, const std::array<float, order__ + 1>& yweights) {
        m_xweights = xweights;
        m_yweights = yweights;
    }

    IIRFilter() = default;

    float push_sample(const float& input) {
        m_input_buffer.push_sample(input);
        m_output_buffer.shift_buffer();

        update_filter();
        return m_output_buffer[0];
    }

    float get_output() const {
        return m_output_buffer[0];
    }

    void set_weights(const std::array<float, order__ + 1>& xweights, const std::array<float, order__ + 1>& yweights) {
        m_xweights = xweights;
        m_yweights = yweights;
    }

    void reset_filter() {
        m_input_buffer.clear_buffer();
        m_output_buffer.clear_buffer();
    }

    float freq_response_magnitude(float fr) {
        float num = zpoly_magnitude<order__>(m_xweights, 2 * std::numbers::pi * fr);
        float den = zpoly_magnitude<order__>(m_yweights, 2 * std::numbers::pi * fr);

        return num / den;
    }

private:
    void update_filter() {
        float x_wsum = mulacc_no_simd<order__ + 1, 0>(m_input_buffer.as_array(), m_xweights); // X * xi where i ranges from 0 to order
        float y_wsum = mulacc_no_simd<order__ + 1, 1>(m_output_buffer.as_array(), m_yweights); // Y * yi where i ranges from 1 to order

        float new_sample = (x_wsum - y_wsum);
        m_output_buffer[0] = new_sample;
    }

    SampleBuffer<order__ + 1> m_input_buffer;
    SampleBuffer<order__ + 1> m_output_buffer;

    std::array<float, order__ + 1> m_xweights;
    std::array<float, order__ + 1> m_yweights;
};


#endif //OPENDSP_IIRFILTER_H
