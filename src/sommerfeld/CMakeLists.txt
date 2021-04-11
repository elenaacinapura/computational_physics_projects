add_executable(sommerfeld sommerfeld.c util.c)

target_compile_options(sommerfeld PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(sommerfeld compphy-lib)

add_custom_target(run-sommerfeld
    COMMAND sommerfeld
    DEPENDS sommerfeld
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)