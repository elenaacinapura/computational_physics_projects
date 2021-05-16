add_executable(bound-numerov-periodic periodic.c ../routines.c)

target_compile_options(bound-numerov-periodic PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(bound-numerov-periodic compphy-lib)

add_custom_target(run-bound-numerov-periodic
    COMMAND bound-numerov-periodic
    DEPENDS bound-numerov-periodic
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)