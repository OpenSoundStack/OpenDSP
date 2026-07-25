#include "precise_acc.h"

PreciseAcc::PreciseAcc() {
    m_acc = 0.0;
    m_compensation = 0.0;
    m_t = 0.0;
    m_y = 0.0;
}

double PreciseAcc::get_accumulator() const {
    return m_acc;
}

void PreciseAcc::operator+=(double v) {
    return accumulate(v);
}

void PreciseAcc::operator-=(double v) {
    return accumulate(v);
}

void PreciseAcc::operator+=(float v) {
    return accumulate(v);
}

void PreciseAcc::operator-=(float v) {
    return accumulate(v);
}

PreciseAcc& PreciseAcc::operator=(float v) {
    set_acc(v);

    return *this;
}

PreciseAcc& PreciseAcc::operator=(double v) {
    set_acc(v);
    return *this;
}

PreciseAcc::operator float() const {
    return static_cast<float>(m_acc);
}

PreciseAcc::operator double() const {
    return m_acc;
}
