#========================= COSH-4 ==============================
add_executable(bound-numerov-cosh cosh.c ../routines.c)

target_compile_options(bound-numerov-cosh PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(bound-numerov-cosh compphy-lib)

add_custom_target(run-bound-numerov-cosh
    COMMAND bound-numerov-cosh
    DEPENDS bound-numerov-cosh
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

#========================= LENNARD JONES ==============================
add_executable(bound-numerov-lj lennard_jones.c ../routines.c)

target_compile_options(bound-numerov-lj PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(bound-numerov-lj compphy-lib)

add_custom_target(run-bound-numerov-lj
    COMMAND bound-numerov-lj
    DEPENDS bound-numerov-lj
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)