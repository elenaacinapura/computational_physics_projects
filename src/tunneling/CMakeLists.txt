# Calculation of the solution to the schrodinger equation
add_executable(solve_equation solve_equation.c util.c)

target_compile_options(solve_equation PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(solve_equation compphy-lib)

add_custom_target(run-tunneling-solve
    COMMAND solve_equation
    DEPENDS solve_equation
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)


# Calculation of T(E) for a square potential
add_executable(transmission transmission.c util.c)

target_compile_options(transmission PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(transmission compphy-lib)

add_custom_target(run-tunneling-transmission
    COMMAND transmission
    DEPENDS transmission
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)