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

#ifndef OPENDSP_SAMPLE_BUFFER_H
#define OPENDSP_SAMPLE_BUFFER_H

#include <cstring>
#include <ranges>

template<int buflen__>
class SampleBuffer {
public:
    SampleBuffer() {
        clear_buffer();
        m_buffer_range = std::ranges::subrange(m_buffer.begin(), m_buffer.end());
    }

    ~SampleBuffer() = default;

    void shift_buffer() {
        auto* data_ptr = m_buffer.data();
        std::memcpy(data_ptr + 1, data_ptr, sizeof(float) * (buflen__ - 1));
    }

    void push_sample(const float& sample) {
        shift_buffer();
        m_buffer[0] = sample;
    }

    void clear_buffer() {
        // Zeroing the buffer
        std::memset(m_buffer.data(), 0, m_buffer.size());
    }

    float operator[](const int idx) const {
        return m_buffer[idx];
    }

    float& operator[](const int idx) {
        return m_buffer[idx];
    }

    std::array<float, buflen__>& as_array() {
        return m_buffer;
    }

    std::ranges::subrange<float*>& as_subrange() {
        return m_buffer_range;
    }

private:
    std::array<float, buflen__> m_buffer;
    std::ranges::subrange<float*> m_buffer_range;
};

#endif //OPENDSP_SAMPLE_BUFFER_H
