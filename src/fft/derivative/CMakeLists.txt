add_executable(derivative derivative.c)

target_compile_options(derivative PRIVATE -O3 -march=native -ffast-math)

target_link_libraries(derivative compphy-lib gsl)

add_custom_target(run-derivative
    COMMAND derivative
    DEPENDS derivative
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/derivative.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(plot-derivative
    COMMAND python derivative.py 2>&1 | grep -v 'libGL'
    DEPENDS derivative.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)