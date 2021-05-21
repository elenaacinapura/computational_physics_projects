add_executable(ljscatt ljscatt.c util.c)

target_compile_options(ljscatt PRIVATE -Wall -O3 -march=native -ffast-math -lgsl)

target_link_libraries(ljscatt m compphy-lib gsl)

add_custom_target(run-ljscatt
    COMMAND ljscatt
    DEPENDS ljscatt 
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_sigma.gp
DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(plot-ljscatt
    COMMAND gnuplot plot_sigma.gp
    DEPENDS plot_sigma.gp
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)