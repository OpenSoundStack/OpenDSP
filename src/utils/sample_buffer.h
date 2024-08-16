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
