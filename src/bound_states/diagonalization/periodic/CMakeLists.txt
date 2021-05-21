add_executable(bound-diag-periodic periodic.c)

target_compile_options(bound-diag-periodic PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(bound-diag-periodic compphy-lib)

add_custom_target(run-bound-diag-periodic
    COMMAND bound-diag-periodic 
    DEPENDS bound-diag-periodic
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/periodic.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(plot-bound-diag-periodic
    COMMAND python periodic.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)