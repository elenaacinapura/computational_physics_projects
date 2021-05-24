#======================= GROUND STATE ===========================
# numerov
add_executable(tunneling-ground-numerov ground_numerov.c routines.c)

target_compile_options(tunneling-ground-numerov PRIVATE -O3 -march=native)

target_link_libraries(tunneling-ground-numerov compphy-lib)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ground_numerov.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)
     
add_custom_target(run-es2-1-numerov
    COMMAND tunneling-ground-numerov && python ground_numerov.py 2>&1 | grep -v 'libGL'
    DEPENDS tunneling-ground-numerov
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
    
# diag
add_executable(tunneling-ground-diag ground_diag.c)

target_compile_options(tunneling-ground-diag PRIVATE -O3 -march=native)

target_link_libraries(tunneling-ground-diag compphy-lib)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/ground_diag.py
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)
add_custom_target(run-es2-1-diag
    COMMAND tunneling-ground-diag && python ground_diag.py 2>&1 | grep -v 'libGL'
    DEPENDS tunneling-ground-diag
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

#======================= EVOLUTIION ===========================

add_executable(tunneling-evol1 evolution.c)

target_compile_options(tunneling-evol1 PRIVATE -O3 -march=native)

target_link_libraries(tunneling-evol1 m compphy-lib gsl openblas)

add_custom_target(run-es2-2
    COMMAND tunneling-evol1
    DEPENDS tunneling-evol1
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_animation.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(plot-es2-2-animation
     COMMAND gnuplot plot_animation.gp
     DEPENDS plot_animation.gp
     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/probability.py
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(plot-es2-2-prob
     COMMAND python probability.py
     DEPENDS probability.py
     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
