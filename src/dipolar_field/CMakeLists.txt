add_executable(dipolar dipolar.c)

target_compile_options(dipolar PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(dipolar compphy-lib)

add_custom_target(run-dipolar
    COMMAND dipolar
    DEPENDS dipolar
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)