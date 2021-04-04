add_executable(tunneling tunneling.c util.c)

target_compile_options(tunneling PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(tunneling compphy-lib)

add_custom_target(run-tunneling
    COMMAND tunneling
    DEPENDS tunneling
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)
