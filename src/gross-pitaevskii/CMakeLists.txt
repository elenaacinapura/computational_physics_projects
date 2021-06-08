add_executable(gross gp.c)

target_compile_options(gross PRIVATE -O3 -march=native)

target_link_libraries(gross compphy-lib gsl)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/gp.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(run-gross-pitaevskii
    COMMAND ulimit -s unlimited && ./gross && python gp.py 2>&1 | grep -v 'libGL'
    DEPENDS gross
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

add_custom_target(plot-gross-pitaevskii
    COMMAND python gp.py 2>&1 | grep -v 'libGL'
    DEPENDS gp.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)