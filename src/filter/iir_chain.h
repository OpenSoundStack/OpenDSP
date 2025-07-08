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

#ifndef IIR_CHAIN_H
#define IIR_CHAIN_H

#include "simd/simd_op.h"
#include "simd/matmath.h"
#include "utils/sample_buffer.h"
#include "iirfilter.h"

template<int order__>
class IIRChain {
public:
    IIRChain() {
        static_assert(order__ <= 4 && order__ > 0);

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                m_xweights.coefs[i][j] = 0;
                m_yweights.coefs[i][j] = 0;
            }
        }

        m_filter_count = 0;
    }

    ~IIRChain() = default;

    float push_sample(const float& input) {
        m_input_buffer.push_sample(input);

        for (auto& ob : m_output_buffers) {
            ob.shift_buffer();
        }

        update_chain();
        return m_chain_out;
    }

    void add_filter(const IIRFilter<order__>& filter) {
        assert(m_filter_count <= 4);

        auto x = filter.get_xweights();
        auto y = filter.get_yweights();

        memcpy(m_xweights.coefs[m_filter_count], x.data(), sizeof(float) * (order__ + 1));
        memcpy(m_yweights.coefs[m_filter_count], y.data(), sizeof(float) * (order__ + 1));

        m_filter_count++;
    }

private:
    void update_chain() {
        Vec4 chain_x = mat4_vec4_mul(m_xweights, m_input_buffer.get_buffer());
        Vec4 chain_y = {0};

        for (int i = 0; i < m_filter_count; i++) {
            chain_y.elems[i] = mulacc<order__ + 1, 1>(m_output_buffers[i].get_buffer(), m_yweights.coefs[i]);
        }

        Vec4 result_vector = vec4_diff(chain_x, chain_y);

        for (int i = 0; i < m_filter_count; i++) {
            m_output_buffers[i][0] = result_vector.elems[i];
        }

        float new_sample = 1.0f;
        for (int i = 0; i < m_filter_count; i++) {
            new_sample *= result_vector.elems[i];
        }

        m_chain_out = new_sample;
        //m_output_buffer[0] = new_sample;
    }

    Mat4x4 m_xweights;
    Mat4x4 m_yweights;

    SampleBuffer<order__ + 1> m_input_buffer;
    SampleBuffer<order__ + 1> m_output_buffers[4];

    int m_filter_count;
    float m_chain_out;
};

#endif //IIR_CHAIN_H
