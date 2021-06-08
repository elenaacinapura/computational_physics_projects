add_executable(schrodinger schrodinger.c)

target_compile_options(schrodinger PRIVATE -O3 -march=native)

target_link_libraries(schrodinger compphy-lib gsl)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/schrodinger.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(run-schrodinger-evolution
    COMMAND schrodinger && python schrodinger.py 2>&1 | grep -v 'libGL'
    DEPENDS schrodinger
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

add_custom_target(plot-schrodinger-evolution
    COMMAND python schrodinger.py 2>&1 | grep -v 'libGL'
    DEPENDS schrodinger.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)