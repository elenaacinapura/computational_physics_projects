add_executable(sir sir.c)

target_compile_options(sir PRIVATE -O3 -march=native)

target_link_libraries(sir compphy-lib gsl)

add_custom_target(run-sir
    COMMAND sir
    DEPENDS sir
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/sir.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(plot-sir
    COMMAND python sir.py 2>&1 | grep -v 'libGL'
    DEPENDS sir.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)