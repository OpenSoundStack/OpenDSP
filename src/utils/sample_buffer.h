// This file is part of the Open Audio Live System project, a live audio environment
// Copyright (c) 2026 - Mathis DELGADO
//
// This project is distributed under the Creative Commons CC-BY-NC-SA licence. https://creativecommons.org/licenses/by-nc-sa/4.0

#ifndef OPENDSP_SAMPLE_BUFFER_H
#define OPENDSP_SAMPLE_BUFFER_H

#include <cstring>
#include <ranges>

template<int buflen__>
class SampleBuffer {
public:
    SampleBuffer() {
        clear_buffer();
    }

    ~SampleBuffer() = default;

    void shift_buffer() {
        std::memcpy(m_buffer + 1, m_buffer, sizeof(float) * (buflen__ - 1));
    }

    void push_sample(const float& sample) {
        shift_buffer();
        m_buffer[0] = sample;
    }

    void clear_buffer() {
        // Zeroing the buffer
        std::memset(m_buffer, 0, sizeof(float) * buflen__);
    }

    float operator[](const int idx) const {
        return m_buffer[idx];
    }

    float& operator[](const int idx) {
        return m_buffer[idx];
    }

    float* get_buffer() {
        return m_buffer;
    }

private:
    float m_buffer[buflen__];
};

#endif //OPENDSP_SAMPLE_BUFFER_H
