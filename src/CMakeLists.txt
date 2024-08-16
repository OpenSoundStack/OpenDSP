add_library(OpenDSP SHARED
        dsp.cpp
        simd/simd_op.h
        filter/iirfilter.h
        utils/sample_buffer.h
        filter/analog/lowpass.h
        utils/values.h
        utils/values.cpp
        filter/analog/highpass.h
)

target_include_directories(OpenDSP PUBLIC ${CMAKE_SOURCE_DIR}/src)