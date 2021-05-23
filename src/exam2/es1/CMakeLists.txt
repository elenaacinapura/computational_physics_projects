add_executable(ljscatt ljscatt.c util.c)

target_compile_options(ljscatt PRIVATE -Wall -O3 -march=native -ffast-math -lgsl)

target_link_libraries(ljscatt m compphy-lib gsl)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scattering.py
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(run-es1
    COMMAND ljscatt && python scattering.py 2>&1 | grep -v 'libGL'
    DEPENDS ljscatt 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
