// This file is part of the Open Audio Live System project, a live audio environment
// Copyright (c) 2026 - Mathis DELGADO
//
// This project is distributed under the Creative Commons CC-BY-NC-SA licence. https://creativecommons.org/licenses/by-nc-sa/4.0

#include "enveloppe.h"

Enveloppe::Enveloppe(int time_constant, int sampling_rate) {
    m_time_constant = time_constant;
    m_sampling_rate = sampling_rate;
    m_delay_samples = 0;
    m_hold_counter = 0;

    init_buffers();
}

Enveloppe::~Enveloppe() {

}

void Enveloppe::set_time_constant(int new_attack) {
    m_time_constant = new_attack;
}

/**
 * @param sample New sample
 * @return Squared mean of the signal
 */
float Enveloppe::push_sample(float sample) {
    update_buffer(m_enveloppe_buffer, sample);
    return  (float)m_accumulator / (float)(m_enveloppe_buffer.size());
}

void Enveloppe::init_buffers() {
    for (int i = 0; i < m_time_constant * (m_sampling_rate / 1000); i++) {
        m_enveloppe_buffer.push(0.0f);
    }

    m_accumulator = 0.0f;
}

void Enveloppe::update_buffer(std::queue<float> &buffer, float new_sample) {
    float sample2 = new_sample * new_sample;

    buffer.push(sample2);
    float first_val = buffer.front();
    buffer.pop();

    double x = sample2 - first_val;
    m_accumulator += x;
}
