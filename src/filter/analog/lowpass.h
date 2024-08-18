#ifndef OPENDSP_LOWPASS_H
#define OPENDSP_LOWPASS_H

#include "filter/iirfilter.h"
#include "utils/values.h"

template<class Tsample>
class LPF_1ord {
public:
    LPF_1ord(const float fc, const float sampling_freq) {
        m_fc = fc;
        m_sampling_freq = sampling_freq;

        m_alpha = compute_alpha(m_fc, m_sampling_freq);

        init_filter();
    }

    ~LPF_1ord() = default;

    Tsample push_sample(const Tsample& s) {
        return m_filter.push_sample(s);
    }

    Tsample get_output() const {
        return m_filter.get_output();
    }

    void set_cutoff(const float fc) {
        m_fc = fc;
        m_alpha = compute_alpha(m_fc, m_sampling_freq);

        update_filter();
    }

    IIRFilter<Tsample, 1>& get_filter() {
        return m_filter;
    }

private:
    void init_filter() {
        auto weights = compute_weights();
        m_filter = IIRFilter<Tsample, 1>(weights[0], weights[1]);
    }

    void update_filter() {
        auto weights = compute_weights();
        m_filter.set_weights(weights[0], weights[1]);
    }

    std::array<std::array<float, 2>, 2> compute_weights() {
        float AB = 1.0f / (1.0f + m_alpha);

        std::array<float, 2> xweights = {
                AB,
                AB
        };

        std::array<float, 2> yweights = {
                1.0f,
                (1 - m_alpha) / (1 + m_alpha)
        };

        return {xweights, yweights};
    }

    float m_alpha;
    float m_fc; // Cutoff frequency
    float m_sampling_freq;

    IIRFilter<Tsample, 1> m_filter;
};

template<class Tsample>
class LPF_2ord {
public:
    LPF_2ord(float fc, float m, float sampling_freq) {
        m_fc = fc;
        m_m = m;
        m_sampling_freq = sampling_freq;

        m_alpha = compute_alpha(m_fc, m_sampling_freq);

        init_filter();
    }

    ~LPF_2ord() = default;

    Tsample push_sample(const Tsample& s) {
        return m_filter.push_sample(s);
    }

    Tsample get_output() const {
        return m_filter.get_output();
    }

    void set_cutoff(float fc) {
        m_fc = fc;
        m_alpha = compute_alpha(m_fc, m_sampling_freq);

        update_filter();
    }

    void set_damping_coef(float m) {
        m_m = m;
        update_filter();
    }

    IIRFilter<Tsample, 2>& get_filter() {
        return m_filter;
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
        float alpha2 = m_alpha * m_alpha;
        float inv_com_den = 1.0 + (2.0 * m_alpha * m_m) + alpha2;
        inv_com_den = 1.0f / inv_com_den;

        std::array<float, 3> xweights = {
                1.0f * inv_com_den,
                2.0f * inv_com_den,
                1.0f * inv_com_den
        };

        std::array<float, 3> yweights = {
                1.0f,
                (2.0f - (2.0f * alpha2)) * inv_com_den,
                (1.0f - (2.0f * m_alpha * m_m) + alpha2) * inv_com_den
        };

        return {xweights, yweights};
    }

    float m_alpha;
    float m_m; // Damping coefficient; Q = 1/(2*m)
    float m_fc; // Cutoff frequency
    float m_sampling_freq;

    IIRFilter<Tsample, 2> m_filter;
};

#endif //OPENDSP_LOWPASS_H
