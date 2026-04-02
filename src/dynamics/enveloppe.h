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

#ifndef OALIVESYSTEM_ENVELOPPE_H
#define OALIVESYSTEM_ENVELOPPE_H

#include <list>
#include <cmath>

class Enveloppe {
public:
    Enveloppe(int time_constant, int sampling_rate);
    ~Enveloppe();

    float push_sample(float sample);

    void set_time_constant(int new_attack);
private:
    void init_buffers();
    void update_buffer(std::list<float>& buffer, float& acc, float new_sample);

    std::list<float> m_enveloppe_buffer;
    float m_enveloppe_acc;

    int m_time_constant;
    int m_sampling_rate;
    int m_delay_samples;

    int m_hold_counter;
};

#endif //OALIVESYSTEM_ENVELOPPE_H