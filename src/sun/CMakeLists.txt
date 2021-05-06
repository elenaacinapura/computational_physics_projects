add_executable(sun sun.c util.c)

target_compile_options(sun PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(sun m compphy-lib)

add_custom_target(run-sun
    COMMAND sun #> out.txt
    DEPENDS sun
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)