#ifndef OPENDSP_HIGHPASS_H
#define OPENDSP_HIGHPASS_H

#include "filter/iirfilter.h"
#include "utils/values.h"

template<class Tsample>
class HPF_1ord {
public:
    HPF_1ord(float fc, float sampling_freq) {
        m_fc = fc;
        m_sampling_freq = sampling_freq;
        m_alpha = compute_alpha(fc, sampling_freq);

        init_filter();
    }
    ~HPF_1ord() = default;

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

    IIRFilter<Tsample, 1>& get_filter() {
        return m_filter;
    }

private:
    void init_filter() {
        auto weights = compute_weights();
        m_filter = IIRFilter<Tsample, 1>{weights[0], weights[1]};
    }

    void update_filter() {
        auto weights = compute_weights();
        m_filter.set_weights(weights[0], weights[1]);
    }

    std::array<std::array<float, 2>, 2> compute_weights() {
        float A = m_alpha / (1.0f + m_alpha);

        std::array<float, 2> xweights = {
                A,
                -A
        };
        std::array<float, 2> yweights = {
                1.0f,
                (1.0f - m_alpha) / (1.0f + m_alpha)
        };

        return {xweights, yweights};
    }

    float m_fc;
    float m_sampling_freq;
    float m_alpha;

    IIRFilter<Tsample, 1> m_filter;
};

template<class Tsample>
class HPF_2ord {
public:
    HPF_2ord(float fc, float m, float sampling_freq) {
        m_fc = fc;
        m_m = m;
        m_sampling_freq = sampling_freq;
        m_alpha = compute_alpha(fc, sampling_freq);

        init_filter();
    }
    ~HPF_2ord() = default;

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
        m_filter = IIRFilter<Tsample, 2>{weights[0], weights[1]};
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
                alpha2 * inv_com_den,
                -2.0f * alpha2 * inv_com_den,
                alpha2 * inv_com_den
        };
        std::array<float, 3> yweights = {
                1.0f,
                (2.0f - 2.0f * alpha2) * inv_com_den,
                (1 - 2 * m_alpha * m_m + alpha2) * inv_com_den
        };

        return {xweights, yweights};
    }

    float m_fc;
    float m_sampling_freq;
    float m_m;
    float m_alpha;

    IIRFilter<Tsample, 2> m_filter;
};

#endif //OPENDSP_HIGHPASS_H
