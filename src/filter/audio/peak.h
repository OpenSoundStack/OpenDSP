#ifndef OPENDSP_PEAK_H
#define OPENDSP_PEAK_H

#include "filter/iirfilter.h"
#include "utils/values.h"

template<class Tsample>
class PeakFilter {
public:
    PeakFilter(float fc, float Q, float gain, float sampling_freq) {
        m_fc = fc;
        m_Q = Q;
        m_gain = gain;
        m_sampling_freq = sampling_freq;

        init_filter();
    }

    ~PeakFilter() = default;

    Tsample push_sample(const Tsample& s) {
        return m_filter.push_sample(s);
    }

    Tsample get_output() const {
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

private:
    void init_filter() {
        auto weights = compute_weights();
        m_filter = IIRFilter<Tsample, 2>(weights[0], weights[1]);
    }

    void update_filter() {
        auto weights = compute_weights();
        m_filter.set_weights(weights[0], weights[1]);
    }

    std::array<std::array<float, 3>, 2> compute_weights() {
        float A = std::pow(10.0f, m_gain / 40.0f);
        float rQ = m_Q / A;
        float w0 = 2.0f * std::numbers::pi * (m_fc / m_sampling_freq);
        float alpha = std::sin(w0) / (2.0f * rQ);

        float inv_com_den = 1.0f + (alpha / A);

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

    IIRFilter<Tsample, 2> m_filter;
};

#endif //OPENDSP_PEAK_H
