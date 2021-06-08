add_executable(fluid fluid.c util.c FCC.c)

target_compile_options(fluid PRIVATE -O3 -march=native -ffast-math)

target_link_libraries(fluid compphy-lib)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/energies.gp
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/T.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/g.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(run-molecular-dynamics
    COMMAND fluid && gnuplot energies.gp && gnuplot T.gp
    DEPENDS fluid
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
