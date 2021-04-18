add_executable(bound bound_states.c util.c)

target_compile_options(bound PRIVATE -O3 -march=native -ffast-math)

target_link_libraries(bound compphy-lib)

add_custom_target(run-bound
    COMMAND bound
    DEPENDS bound
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/delta.gp
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})