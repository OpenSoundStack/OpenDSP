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
    : m_signal_enveloppe_attack(attack_ms, sampling_rate),
    m_signal_enveloppe_release(release_ms, sampling_rate)
{
    m_attack_ms = attack_ms;
    m_release_ms = release_ms;
    m_hold_ms = hold_ms;

    m_last_env_value = 0.0f;
    m_current_env_value = 0.0f;
    m_deriv_counter = 0;
    m_last_attack_value = 1.0f;
    m_last_release_value = 1.0f;

    m_hold_time_sample = m_hold_ms * (sampling_rate / 1000);
    m_hold_counter = 0;

    m_transfer_function = std::move(transfer_function);

    m_state = DynamicState::DYN_RELEASE;

    init_delay_buffer();
}

Dynamics::~Dynamics() {

}

float Dynamics::push_sample(float sample) {
    float sq_mean_attack = m_signal_enveloppe_attack.push_sample(sample);
    float level_lin_attack = std::sqrt(sq_mean_attack);

    float sq_mean_release = m_signal_enveloppe_release.push_sample(sample);
    float level_lin_release = process_delay(std::sqrt(sq_mean_release));

    float reduction_attack = 1.0f;
    if (level_lin_attack != 0.0f) {
        reduction_attack = m_transfer_function(level_lin_attack) / level_lin_attack;
    }

    float reduction_release = 1.0f;
    if (level_lin_release != 0.0f) {
        reduction_release = m_transfer_function(level_lin_release) / level_lin_release;
    }

    float diff_env = differentiate_enveloppe(level_lin_attack);
    constexpr float hysteresis = 0.001f;

    switch (m_state) {
        case DynamicState::DYN_ATTACK:
            if (diff_env < -hysteresis) {
                m_state = DynamicState::DYN_HOLD;
                m_hold_counter = m_hold_time_sample;
            }
            break;
        case DynamicState::DYN_RELEASE:
            if (diff_env > hysteresis) {
                m_state = DynamicState::DYN_ATTACK;
            }
            break;
        case DynamicState::DYN_HOLD:
            if (m_hold_counter == 0) {
                m_state = DynamicState::DYN_RELEASE;
            }
            break;
        default:
            break;
    }

    float reduction = 1.0f;

    switch (m_state) {
        case DynamicState::DYN_ATTACK:
            reduction = reduction_attack;
            m_last_attack_value = reduction_attack;
            break;
        case DynamicState::DYN_RELEASE:
            reduction = reduction_release;
            break;
        case DynamicState::DYN_LOCKED:
            reduction = m_last_attack_value;
            break;
        case DynamicState::DYN_HOLD:
            reduction = m_last_attack_value;
            m_hold_counter--;
            break;
        default:
            break;
    }

    return reduction;
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
