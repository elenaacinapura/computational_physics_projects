# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /media/Dati/Git/ComputationalPhysics

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /media/Dati/Git/ComputationalPhysics/extern

# Utility rule file for run-fluid.

# Include the progress variables for this target.
include src/lennard_jones_fluid/CMakeFiles/run-fluid.dir/progress.make

src/lennard_jones_fluid/CMakeFiles/run-fluid: src/lennard_jones_fluid/fluid
	cd /media/Dati/Git/ComputationalPhysics/extern/src/lennard_jones_fluid && ./fluid

run-fluid: src/lennard_jones_fluid/CMakeFiles/run-fluid
run-fluid: src/lennard_jones_fluid/CMakeFiles/run-fluid.dir/build.make

.PHONY : run-fluid

# Rule to build all files generated by this target.
src/lennard_jones_fluid/CMakeFiles/run-fluid.dir/build: run-fluid

.PHONY : src/lennard_jones_fluid/CMakeFiles/run-fluid.dir/build

src/lennard_jones_fluid/CMakeFiles/run-fluid.dir/clean:
	cd /media/Dati/Git/ComputationalPhysics/extern/src/lennard_jones_fluid && $(CMAKE_COMMAND) -P CMakeFiles/run-fluid.dir/cmake_clean.cmake
.PHONY : src/lennard_jones_fluid/CMakeFiles/run-fluid.dir/clean

src/lennard_jones_fluid/CMakeFiles/run-fluid.dir/depend:
	cd /media/Dati/Git/ComputationalPhysics/extern && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /media/Dati/Git/ComputationalPhysics /media/Dati/Git/ComputationalPhysics/src/lennard_jones_fluid /media/Dati/Git/ComputationalPhysics/extern /media/Dati/Git/ComputationalPhysics/extern/src/lennard_jones_fluid /media/Dati/Git/ComputationalPhysics/extern/src/lennard_jones_fluid/CMakeFiles/run-fluid.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/lennard_jones_fluid/CMakeFiles/run-fluid.dir/depend

