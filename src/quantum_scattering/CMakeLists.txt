add_executable(quantum_scattering quantum_scattering.c util.c)

target_compile_options(quantum_scattering PRIVATE -Wall -O3 -march=native -ffast-math -lgsl)

target_link_libraries(quantum_scattering compphy-lib gsl)

add_custom_target(run-quantum-scattering
    COMMAND quantum_scattering
    DEPENDS quantum_scattering
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/sigma_tot.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})