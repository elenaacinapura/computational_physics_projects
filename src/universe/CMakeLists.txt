add_executable(universe universe.c util.c)

target_compile_options(universe PRIVATE -Wall -O3 -march=native -ffast-math)
# target_link_options(universe PRIVATE -fsanitize=address)

target_link_libraries(universe compphy-lib)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/universe.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(run-universe
    COMMAND universe && python universe.py 2> grep -v 'libGL' 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
