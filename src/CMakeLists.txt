add_library(OpenDSP SHARED
        dsp.cpp
        simd/simd_op.h
        filter/iirfilter.h
        utils/sample_buffer.h
        filter/analog/lowpass.h
        utils/values.h
        utils/values.cpp
        filter/analog/highpass.h
        filter/analog/bandpass.h
        filter/audio/peak.h
        simd/matmath.h
        filter/iir_chain.h
        dynamics/enveloppe.h
        dynamics/enveloppe.cpp
        dynamics/dynamics.h
        dynamics/dynamics.cpp
)

target_include_directories(OpenDSP PUBLIC ${PROJECT_SOURCE_DIR}/src)