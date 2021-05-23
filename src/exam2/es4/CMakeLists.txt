add_executable(es4 percus_yevick.c)

target_compile_options(es4 PRIVATE -O3 -march=native)

target_link_libraries(es4 compphy-lib gsl)


file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/percus_yevick.py
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(run-es4
    COMMAND es4 && python percus_yevick.py 2>&1 | grep -v 'libGL'
    DEPENDS es4
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
