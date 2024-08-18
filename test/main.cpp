#include <iostream>
#include <chrono>
#include <cmath>

#include "filter/analog/lowpass.h"
#include "filter/analog/highpass.h"
#include "filter/analog/bandpass.h"
#include "filter/audio/peak.h"

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

    BPF_2ord<float> bpf{ 50.0f, 5.0f, 48000.0f };

    PeakFilter<float> peak{ 5000.0f, 10.0f, 5.0f, 48000.0f };
    PeakFilter<float> peak2{ 200.0f, 0.3f, 3.0f, 48000.0f };

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

    int mag_npoints = 150;

    std::vector<float> frmag{};
    frmag.reserve(mag_npoints);

    auto start = high_resolution_clock::now();

    for(int i = 0; i < mag_npoints; i++) {
        float fr = i / (2.0f * mag_npoints);

        float mag1 = peak.get_filter().freq_response_magnitude(fr);
        float mag2 = peak2.get_filter().freq_response_magnitude(fr);

        float dbMag = 20.0f * log10(mag1 * mag2);

        frmag.push_back(dbMag);
    }

    auto stop = high_resolution_clock::now();
    auto dur = stop - start;
    auto viz_gen = duration_cast<microseconds>(dur).count();

    std::cout << "FRMAG duration : " << viz_gen << " us" << std::endl;

    for(auto e : sigin) {
        sigout.push_back(peak.push_sample(e));
    }

    plt::plot(frmag);
    plt::show();

    plt::plot(sigout);
    plt::show();

    return 0;
}