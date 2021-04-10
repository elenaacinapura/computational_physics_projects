cmake_minimum_required(VERSION 3.13)
project(ComputationalPhysicsProject)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_subdirectory(extern/computational_physics_lib)
# List here the directories in src
add_subdirectory(src/lennard_jones_fluid)
add_subdirectory(src/tunneling)
add_subdirectory(src/sommerfeld)
add_subdirectory(src/variational_method)
add_subdirectory(src/dipolar_field)