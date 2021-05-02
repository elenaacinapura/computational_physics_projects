add_executable(gravity gravity.c)

target_compile_options(gravity PRIVATE -O3 -march=native -ffast-math)

target_link_libraries(gravity compphy-lib gsl)

add_custom_target(run-gravity
    COMMAND gravity
    DEPENDS gravity
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/gravity.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_target(plot-gravity
    COMMAND python gravity.py
    DEPENDS gravity.py
    WORKING_DIRECTORY ${CMAKE_PROJECT_DIR}
)