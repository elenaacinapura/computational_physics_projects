add_executable(variational variational.c util.c)

target_compile_options(variational PRIVATE -Wall -O3 -march=native -ffast-math)

target_link_libraries(variational compphy-lib)

add_custom_target(run-variational
    COMMAND variational
    DEPENDS variational
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_E.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})