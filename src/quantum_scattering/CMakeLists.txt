add_executable(quantum_scattering quantum_scattering.c util.c)

target_compile_options(quantum_scattering PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(quantum_scattering compphy-lib)

add_custom_target(run-quantum-scattering
    COMMAND quantum_scattering
    DEPENDS quantum_scattering
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)