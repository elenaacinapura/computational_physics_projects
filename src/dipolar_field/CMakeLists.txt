add_executable(dipolar dipolar.c util.c)

target_compile_options(dipolar PRIVATE -O3 -march=native -ffast-math)

target_link_libraries(dipolar compphy-lib)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_pi.gp
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_pi_45.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(run-dipolar
    COMMAND dipolar && gnuplot plot_pi.gp && gnuplot plot_pi_45.gp
    DEPENDS dipolar
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

