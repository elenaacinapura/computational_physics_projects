add_executable(gross gp.c)

target_compile_options(gross PRIVATE -O3 -march=native)

target_link_libraries(gross compphy-lib gsl)

add_custom_target(run-gross
    COMMAND ulimit -s unlimited && ./gross
    DEPENDS gross
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/gp.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(plot-gross
    COMMAND python gp.py 2>&1 | grep -v 'libGL'
    DEPENDS gp.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)