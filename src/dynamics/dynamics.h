// This file is part of the Open Audio Live System project, a live audio environment
// Copyright (c) 2026 - Mathis DELGADO
//
// This project is distributed under the Creative Commons CC-BY-NC-SA licence. https://creativecommons.org/licenses/by-nc-sa/4.0

#ifndef OALIVESYSTEM_DYNAMICS_H
#define OALIVESYSTEM_DYNAMICS_H

#include <functional>
#include <list>

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

    void set_attack(int attack_ms);
    void set_release(int release_ms);
    void set_hold(int hold_ms);

private:
    float differentiate_enveloppe(float enveloppe_sample);
    void init_delay_buffer();
    void make_adsr_coefs();

    float process_delay(float sample);
    float adsr_process(float sample);

    void adjust_delay_buffer();

    std::function<float(float)> m_transfer_function;

    Enveloppe m_signal_enveloppe;
    std::list<float> m_delay_buffer;

    int m_attack_ms;
    int m_release_ms;
    int m_hold_ms;

    float m_adsr_att_coef;
    float m_adsr_rel_coef;

    int m_hold_time_sample;
    int m_hold_counter;

    float m_last_env_value;
    float m_current_env_value;
    int m_deriv_counter;

    float m_enveloppe;
    float m_delayed_enveloppe;

    float m_last_attack_value;
    DynamicState m_state;

    int m_sampling_rate;
};

#endif //OALIVESYSTEM_DYNAMICS_H