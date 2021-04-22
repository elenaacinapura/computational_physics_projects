add_executable(bound-delta delta_method.c util.c)

target_compile_options(bound-delta PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(bound-delta compphy-lib)

add_custom_target(run-bound-delta
    COMMAND bound-delta
    DEPENDS bound-delta
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/delta.gp
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

     add_custom_target(plot-bound-delta
     COMMAND gnuplot delta.gp
     DEPENDS delta.gp
     WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

add_executable(bound-diag diag_method.c util.c)

target_compile_options(bound-diag PRIVATE -O3 -march=native -ffast-math)

target_link_libraries(bound-diag compphy-lib)

add_custom_target(run-bound-diag
    COMMAND bound-diag
    DEPENDS bound-diag
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)