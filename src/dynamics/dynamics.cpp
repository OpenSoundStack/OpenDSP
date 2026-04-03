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

#include "dynamics.h"

Dynamics::Dynamics(
    std::function<float(float)> transfer_function,
    int attack_ms, int release_ms, int hold_ms,
    int sampling_rate
)
    : m_signal_enveloppe(10, sampling_rate)
{
    m_attack_ms = attack_ms;
    m_release_ms = release_ms;
    m_hold_ms = hold_ms;

    m_last_env_value = 0.0f;
    m_current_env_value = 0.0f;
    m_deriv_counter = 0;
    m_last_attack_value = 1.0f;
    m_last_release_value = 1.0f;

    m_current_enveloppe = 1.0f;
    m_enveloppe = 0.0f;

    m_hold_time_sample = m_hold_ms * (sampling_rate / 1000);
    m_hold_counter = 0;

    m_transfer_function = std::move(transfer_function);

    m_state = DynamicState::DYN_RELEASE;

    init_delay_buffer();
    make_adsr_coefs();
}

Dynamics::~Dynamics() {

}

void Dynamics::make_adsr_coefs() {
    m_adsr_att_coef = exp(-1.0f / (m_attack_ms * 96.0f));
    m_adsr_rel_coef = exp(-1.0f / (m_release_ms * 96.0f));
}

float Dynamics::adsr_process(float sample) {
    float sq_mean_enveloppe = m_signal_enveloppe.push_sample(sample);
    float level_lin_enveloppe = std::sqrt(sq_mean_enveloppe);

    float exp_coef = m_state == DynamicState::DYN_ATTACK ? m_adsr_att_coef : m_adsr_rel_coef;

    float coef = exp_coef * (m_enveloppe - level_lin_enveloppe);
    m_enveloppe = level_lin_enveloppe + coef;
    m_delayed_enveloppe = process_delay(m_enveloppe);

    return level_lin_enveloppe;
}

float Dynamics::push_sample(float sample) {
    float level_lin_enveloppe = adsr_process(sample);

    float selected_enveloppe = m_state == DynamicState::DYN_RELEASE ? m_delayed_enveloppe : m_enveloppe;

    float transfer_ratio = 1.0f;
    if (level_lin_enveloppe != 0.0f) {
        transfer_ratio = m_transfer_function(selected_enveloppe) / selected_enveloppe;
    }

    constexpr float hysteresis = 0.001f;

    static int time = 0;
    time++;

    float gain = 1.0f;
    switch (m_state) {
        case DynamicState::DYN_ATTACK:
            gain = transfer_ratio;
            break;
        case DynamicState::DYN_RELEASE:
            gain = transfer_ratio;
            break;
        case DynamicState::DYN_HOLD:
            m_hold_counter--;
            gain = m_last_attack_value;

            break;
        default:
            break;
    }

    float diff_env = differentiate_enveloppe(m_enveloppe);
    switch (m_state) {
        case DynamicState::DYN_ATTACK:
            if (diff_env < -hysteresis) {
                m_state = m_hold_ms == 0 ? DynamicState::DYN_RELEASE : DynamicState::DYN_HOLD;
                m_hold_counter = m_hold_time_sample;

                TRACE_STATE(HOLD, time);
            }
            break;
        case DynamicState::DYN_RELEASE:
            if (diff_env > hysteresis) {
                m_state = DynamicState::DYN_ATTACK;
                TRACE_STATE(ATTACK, time);
            }
            break;
        case DynamicState::DYN_HOLD:
            if (m_hold_counter == 0) {
                m_state = DynamicState::DYN_RELEASE;
                TRACE_STATE(RELEASE, time);
            }

            if (diff_env > hysteresis) {
                m_state = DynamicState::DYN_ATTACK;
                TRACE_STATE(ATTACK, time);
            }

            break;
        default:
            break;
    }

    m_last_attack_value = gain;
    return gain;
}

float Dynamics::differentiate_enveloppe(float enveloppe_sample) {
    if (m_deriv_counter == 0) {
        m_last_env_value = m_current_env_value;
        m_current_env_value = enveloppe_sample;
    }

    m_deriv_counter = (m_deriv_counter + 1) % 100;
    return m_current_env_value - m_last_env_value;
}

void Dynamics::init_delay_buffer() {
    m_delay_buffer = std::list<float>(m_hold_time_sample, 0.0f);
}

float Dynamics::process_delay(float sample) {
    m_delay_buffer.push_back(sample);
    float oldest_samples = m_delay_buffer.front();
    m_delay_buffer.pop_front();

    return oldest_samples;
}
