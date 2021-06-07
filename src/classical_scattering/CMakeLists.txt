#=========================== THETA(B) ======================
add_executable(theta theta.c)

target_compile_options(theta PRIVATE -O3 -march=native)

target_link_libraries(theta compphy-lib gsl)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/theta.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(run-theta
    COMMAND theta && python theta.py 2>&1 | grep -v 'libGL'
    DEPENDS theta
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

add_custom_target(plot-theta
    COMMAND python theta.py 2>&1 | grep -v 'libGL'
    DEPENDS theta.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)