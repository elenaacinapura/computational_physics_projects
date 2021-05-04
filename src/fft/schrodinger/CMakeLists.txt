add_executable(schrodinger schrodinger.c)

target_compile_options(schrodinger PRIVATE -O3 -march=native -ffast-math)

target_link_libraries(schrodinger compphy-lib gsl)

add_custom_target(fft-schrodinger
    COMMAND schrodinger
    DEPENDS schrodinger
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/schrodinger.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(plot-fft-schrodinger
    COMMAND python schrodinger.py 2>&1 | grep -v 'libGL'
    DEPENDS schrodinger.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)