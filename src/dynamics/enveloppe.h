// This file is part of the Open Audio Live System project, a live audio environment
// Copyright (c) 2026 - Mathis DELGADO
//
// This project is distributed under the Creative Commons CC-BY-NC-SA licence. https://creativecommons.org/licenses/by-nc-sa/4.0

#ifndef OALIVESYSTEM_ENVELOPPE_H
#define OALIVESYSTEM_ENVELOPPE_H

#include <queue>
#include <cmath>

class Enveloppe {
public:
    Enveloppe(int time_constant, int sampling_rate);
    ~Enveloppe();

    float push_sample(float sample);

    void set_time_constant(int new_attack);
private:
    void init_buffers();
    void update_buffer(std::queue<float>& buffer, float new_sample);

    std::queue<float> m_enveloppe_buffer;

    double m_enveloppe_acc;
    double m_y;
    double m_t;
    double m_compensation;

    int m_time_constant;
    int m_sampling_rate;
    int m_delay_samples;

    int m_hold_counter;
};

#endif //OALIVESYSTEM_ENVELOPPE_H