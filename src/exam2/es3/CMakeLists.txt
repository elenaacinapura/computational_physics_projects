add_executable(es3 sir.c)

target_compile_options(es3 PRIVATE -O3 -march=native)

target_link_libraries(es3 compphy-lib gsl)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/sir.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(run-es3
    COMMAND es3 && python sir.py 2>&1 | grep -v 'libGL'
    DEPENDS es3
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
