add_executable(structure structure_factor.c)

target_compile_options(structure PRIVATE -O3 -march=native -ffast-math)

target_link_libraries(structure compphy-lib gsl)

add_custom_target(run-structure
    COMMAND structure
    DEPENDS structure
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/structure.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(plot-structure
    COMMAND python structure.py 2>&1 | grep -v 'libGL'
    DEPENDS structure.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)