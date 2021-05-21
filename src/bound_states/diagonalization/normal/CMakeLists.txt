#========================= COSH-4 ==============================
add_executable(bound-diag-cosh cosh.c)

target_compile_options(bound-diag-cosh PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(bound-diag-cosh compphy-lib)

add_custom_target(run-bound-diag-cosh
    COMMAND bound-diag-cosh
    DEPENDS bound-diag-cosh
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

#========================= LENNARD-JONES ==============================
add_executable(bound-diag-lj lennard_jones.c)

target_compile_options(bound-diag-lj PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(bound-diag-lj compphy-lib)

add_custom_target(run-bound-diag-lj
    COMMAND bound-diag-lj
    DEPENDS bound-diag-lj
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)