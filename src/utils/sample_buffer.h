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

template<class T, int buflen__>
class SampleBuffer {
public:
    SampleBuffer() {
        clear_buffer();
    }

    ~SampleBuffer() = default;

    void shift_buffer() {
        std::memcpy(m_buffer + 1, m_buffer, sizeof(T) * (buflen__ - 1));
    }

    void push_sample(const T& sample) {
        shift_buffer();
        m_buffer[0] = sample;
    }

    void clear_buffer() {
        // Zeroing the buffer
        std::memset(m_buffer, 0, sizeof(m_buffer));
    }

    T operator[](const int idx) const {
        return m_buffer[idx];
    }

    T& operator[](const int idx) {
        return m_buffer[idx];
    }

private:
    T m_buffer[buflen__];
};

#endif //OPENDSP_SAMPLE_BUFFER_H
