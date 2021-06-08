add_executable(derivative derivative.c)

target_compile_options(derivative PRIVATE -O3 -march=native -ffast-math)

target_link_libraries(derivative compphy-lib gsl)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/derivative.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(run-derivative-fft
    COMMAND derivative && python derivative.py 2>&1 | grep -v 'libGL'
    DEPENDS derivative
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

add_custom_target(plot-derivative-fft
    COMMAND python derivative.py 2>&1 | grep -v 'libGL'
    DEPENDS derivative.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)