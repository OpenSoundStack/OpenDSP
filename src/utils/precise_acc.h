// This file is part of the Open Audio Live System project, a live audio environment
// Copyright (c) 2026 - Mathis DELGADO
//
// This project is distributed under the Creative Commons CC-BY-NC-SA licence. https://creativecommons.org/licenses/by-nc-sa/4.0

#ifndef PRECISEACC_H
#define PRECISEACC_H

class PreciseAcc {
public:
    PreciseAcc();
    ~PreciseAcc() = default;

    void operator+=(double v);
    void operator-=(double v);

    void operator+=(float v);
    void operator-=(float v);

    PreciseAcc& operator=(float v);
    PreciseAcc& operator=(double v);

    operator float() const;
    operator double() const;

    double get_accumulator() const;

private:
    template<typename T>
    void accumulate(T val) {
        m_y = static_cast<double>(val) - m_compensation;
        m_t = m_acc + m_y;
        m_compensation = (m_t - m_acc) - m_y;
        m_acc = m_t;
    }

    template<typename T>
    void set_acc(T val) {
        m_y = 0.0;
        m_t = 0.0;
        m_compensation = 0.0;

        m_acc = static_cast<double>(val);
    }

    double m_acc;
    double m_t;
    double m_y;
    double m_compensation;
};

#endif //PRECISEACC_H
