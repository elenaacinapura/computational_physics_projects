add_executable(fluid fluid.c util.c FCC.c)

target_compile_options(fluid PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(fluid compphy-lib)

add_custom_target(run-fluid
    COMMAND fluid
    DEPENDS fluid
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/energies.gp
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/T.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/g.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})