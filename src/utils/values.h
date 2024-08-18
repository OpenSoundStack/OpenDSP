#ifndef OPENDSP_VALUES_H
#define OPENDSP_VALUES_H

#include <cmath>
#include <numbers>
#include <array>

float compute_alpha(float cutoff, float samplig_freq);

template<int order__>
float zpoly_magnitude(std::array<float, order__ + 1> poly, float w) {
    // Every z-poly in structured in this way :
    // a + b * z^-1 + c * z^-2 + ...
    // Given that z = e^(j*w) we can compute Re(poly) and Im(poly) then compute the magnitude with
    // SQRT(Re(poly)^2 + Im(poly)^2)

    float Re = 0.0f;
    float Im = 0.0f;

    for(int i = 0; i < (order__ + 1); i++) {
        Re += poly[i] * std::cos(i * w);
        Im += poly[i] * std::sin(i * w);
    }

    return std::sqrt(
            (Re * Re) +
            (Im * Im)
    );
}

#endif //OPENDSP_VALUES_H
