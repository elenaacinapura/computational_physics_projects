add_executable(virial virial.c)

target_compile_options(virial PRIVATE -O3 -march=native)

target_link_libraries(virial compphy-lib gsl)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/virial.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target(run-virial
    COMMAND virial && python virial.py 2>&1 | grep -v 'libGL'
    DEPENDS virial
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)