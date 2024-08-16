#include <iostream>
#include <chrono>
#include <cmath>

#include "filter/analog/lowpass.h"
#include "filter/analog/highpass.h"
#include "filter/analog/bandpass.h"

#define WITHOUT_NUMPY
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

using namespace std::chrono;
using namespace std::chrono_literals;

int main() {
    LPF_2ord<float> lpf{ 500.0, 0.707, 48000.0 };
    LPF_1ord<float> lpf1{500.0, 48000.0};

    HPF_1ord<float> hpf{ 5000.0f, 48000.0f };
    HPF_2ord<float> hpf1{ 500.0f, 0.707f, 48000.0f };

    BPF_2ord<float> bpf{ 5000.0f, 10.0f, 48000.0f };

    int npoints = 10000;

    std::vector<float> sigin{};
    sigin.reserve(npoints);

    float f1 = 50;
    float f2 = 1000;
    float f3 = 5000;
    for(int n = 0; n < npoints; n++) {
        sigin.push_back(
                sin(2.0f * 3.141592f * f1 * n * (1/48000.0f)) +
                sin(2.0f * 3.141592f * f2 * n * (1/48000.0f)) +
                sin(2.0f * 3.141592f * f3 * n * (1/48000.0f))
        );
    }

    std::vector<float> sigout{};
    sigout.reserve(npoints);

    auto start = high_resolution_clock::now();

    for(auto e : sigin) {
        sigout.push_back(bpf.push_sample(e));
    }

    //plt::plot(sigin);
    //plt::show();

    plt::plot(sigout);
    plt::show();

    return 0;
}