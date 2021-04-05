# Calculation of the solution to the schrodinger equation
add_executable(tunneling_solutions find_solutions.c util.c)

target_compile_options(tunneling_solutions PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(tunneling_solutions compphy-lib)

add_custom_target(run-tunneling-sol
    COMMAND tunneling_solutions
    DEPENDS tunneling_solutions
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)


# Calculation of T(E) for a square potential
add_executable(tunneling_T_square T_square.c util.c)

target_compile_options(tunneling_T_square PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(tunneling_T_square compphy-lib)

add_custom_target(run-tunneling-square
    COMMAND tunneling_T_square
    DEPENDS tunneling_T_square
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)