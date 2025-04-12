find_package(Python3 COMPONENTS Development)

add_executable(dsp_test
        main.cpp
        matplotlibcpp.h
        ../../engine/main.cpp
        ../../coreui/main.cpp
        ../../coreui/ui/SignalsWindow.cpp
        ../../coreui/ui/SignalsWindow.h
)

target_link_libraries(dsp_test OpenDSP)
target_link_libraries(dsp_test Python3::Python)
target_include_directories(dsp_test PUBLIC ${CMAKE_SOURCE_DIR}/src)