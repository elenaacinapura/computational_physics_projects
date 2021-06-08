add_executable(structure structure_factor.c)

target_compile_options(structure PRIVATE -O3 -march=native -ffast-math)

target_link_libraries(structure compphy-lib gsl)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/structure.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(run-structure-factor
    COMMAND structure && python structure.py 2>&1 | grep -v 'libGL'
    DEPENDS structure
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

add_custom_target(plot-structure-factor
    COMMAND python structure.py 2>&1 | grep -v 'libGL'
    DEPENDS structure.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)