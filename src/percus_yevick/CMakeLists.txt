add_executable(percus percus_yevick.c)

target_compile_options(percus PRIVATE -O3 -march=native -ffast-math)

target_link_libraries(percus compphy-lib gsl)

add_custom_target(run-percus
    COMMAND percus
    DEPENDS percus
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/percus_yevick.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(plot-percus
    COMMAND python percus_yevick.py 2>&1 | grep -v 'libGL'
    DEPENDS percus_yevick.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)