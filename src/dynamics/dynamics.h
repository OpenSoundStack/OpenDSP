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

#ifndef OALIVESYSTEM_DYNAMICS_H
#define OALIVESYSTEM_DYNAMICS_H

#include <functional>
#include <list>
#include <iostream>

#include "enveloppe.h"

enum class DynamicState {
    DYN_RELEASE,
    DYN_ATTACK,
    DYN_LOCKED,
    DYN_HOLD
};

class Dynamics {
public:
    Dynamics(
        std::function<float(float)> transfer_function,
        int attack_ms, int release_ms, int hold_ms,
        int sampling_rate
    );

    ~Dynamics();

     /**
     * Push sample in the system
     * @param sample Signal sample
     * @return Linear gain reduction
     */
    float push_sample(float sample);

private:
    float differentiate_enveloppe(float enveloppe_sample);
    void init_delay_buffer();
    float process_delay(float sample);

    std::function<float(float)> m_transfer_function;

    Enveloppe m_signal_enveloppe_attack;
    Enveloppe m_signal_enveloppe_release;

    std::list<float> m_delay_buffer;

    int m_attack_ms;
    int m_release_ms;
    int m_hold_ms;

    int m_hold_time_sample;
    int m_hold_counter;

    float m_last_env_value;
    float m_current_env_value;
    int m_deriv_counter;

    float m_last_attack_value;
    float m_last_release_value;

    DynamicState m_state;
};

#endif //OALIVESYSTEM_DYNAMICS_H