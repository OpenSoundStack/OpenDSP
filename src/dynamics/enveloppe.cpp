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

#include "enveloppe.h"

Enveloppe::Enveloppe(int time_constant, int sampling_rate) {
    m_time_constant = time_constant;
    m_sampling_rate = sampling_rate;
    m_delay_samples = 0;

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
    update_buffer(m_enveloppe_buffer, m_enveloppe_acc, sample);
    return  m_enveloppe_acc / (float)(m_enveloppe_buffer.size());
}

void Enveloppe::init_buffers() {
    for (int i = 0; i < m_time_constant * (m_sampling_rate / 1000); i++) {
        m_enveloppe_buffer.push_back(0.0f);
    }

    m_enveloppe_acc = 0.0f;
}

void Enveloppe::update_buffer(std::list<float> &buffer, float& acc, float new_sample) {
    float sample2 = new_sample * new_sample;

    buffer.push_back(sample2);
    float first_val = buffer.front();
    buffer.pop_front();

    acc += sample2;
    acc -= first_val;
}
