#ifndef OPENDSP_IIRFILTER_H
#define OPENDSP_IIRFILTER_H

#include <array>

#include "utils/sample_buffer.h"
#include "utils/values.h"

template<class Tsample, int order__>
class IIRFilter {
public:
    IIRFilter(const std::array<float, order__ + 1>& xweights, const std::array<float, order__ + 1>& yweights) {
        m_xweights = xweights;
        m_yweights = yweights;
    }

    IIRFilter() = default;

    Tsample push_sample(const Tsample& input) {
        m_input_buffer.push_sample(input);
        m_output_buffer.shift_buffer();

        update_filter();
        return m_output_buffer[0];
    }

    Tsample get_output() const {
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
        Tsample x_wsum = wsum<order__ + 1, 0>(m_input_buffer, m_xweights); // X * xi where i ranges from 0 to order
        Tsample y_wsum = wsum<order__ + 1, 1>(m_output_buffer, m_yweights); // Y * yi where i ranges from 1 to order

        Tsample new_sample = (x_wsum - y_wsum);
        m_output_buffer[0] = new_sample;
    }

    template<int size__, int init__>
    Tsample wsum(SampleBuffer<Tsample, size__> samples, std::array<float, size__> weights) {
        Tsample accumulator = 0;

        // Samples are ordered in this way : samples[i] = X(n - i)
        // Weights are ordered in the same way : weights[i] = a_i
        // We are doing the sum of X(n - i) * a_i
        for(int i = init__; i < size__; i++) {
            accumulator += samples[i] * weights[i];
        }

        return accumulator;
    }

    SampleBuffer<Tsample, order__ + 1> m_input_buffer;
    SampleBuffer<Tsample, order__ + 1> m_output_buffer;

    std::array<float, order__ + 1> m_xweights;
    std::array<float, order__ + 1> m_yweights;
};


#endif //OPENDSP_IIRFILTER_H
