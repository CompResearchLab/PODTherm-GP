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
CMAKE_SOURCE_DIR = /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/build

# Include any dependencies generated for this target.
include CMakeFiles/calculate_A_Matrix.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/calculate_A_Matrix.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/calculate_A_Matrix.dir/flags.make

CMakeFiles/calculate_A_Matrix.dir/src/main.cpp.o: CMakeFiles/calculate_A_Matrix.dir/flags.make
CMakeFiles/calculate_A_Matrix.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/calculate_A_Matrix.dir/src/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/calculate_A_Matrix.dir/src/main.cpp.o -c /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/src/main.cpp

CMakeFiles/calculate_A_Matrix.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/calculate_A_Matrix.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/src/main.cpp > CMakeFiles/calculate_A_Matrix.dir/src/main.cpp.i

CMakeFiles/calculate_A_Matrix.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/calculate_A_Matrix.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/src/main.cpp -o CMakeFiles/calculate_A_Matrix.dir/src/main.cpp.s

CMakeFiles/calculate_A_Matrix.dir/src/helpers.cpp.o: CMakeFiles/calculate_A_Matrix.dir/flags.make
CMakeFiles/calculate_A_Matrix.dir/src/helpers.cpp.o: ../src/helpers.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/calculate_A_Matrix.dir/src/helpers.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/calculate_A_Matrix.dir/src/helpers.cpp.o -c /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/src/helpers.cpp

CMakeFiles/calculate_A_Matrix.dir/src/helpers.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/calculate_A_Matrix.dir/src/helpers.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/src/helpers.cpp > CMakeFiles/calculate_A_Matrix.dir/src/helpers.cpp.i

CMakeFiles/calculate_A_Matrix.dir/src/helpers.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/calculate_A_Matrix.dir/src/helpers.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/src/helpers.cpp -o CMakeFiles/calculate_A_Matrix.dir/src/helpers.cpp.s

CMakeFiles/calculate_A_Matrix.dir/src/Matrix.cpp.o: CMakeFiles/calculate_A_Matrix.dir/flags.make
CMakeFiles/calculate_A_Matrix.dir/src/Matrix.cpp.o: ../src/Matrix.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/calculate_A_Matrix.dir/src/Matrix.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/calculate_A_Matrix.dir/src/Matrix.cpp.o -c /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/src/Matrix.cpp

CMakeFiles/calculate_A_Matrix.dir/src/Matrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/calculate_A_Matrix.dir/src/Matrix.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/src/Matrix.cpp > CMakeFiles/calculate_A_Matrix.dir/src/Matrix.cpp.i

CMakeFiles/calculate_A_Matrix.dir/src/Matrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/calculate_A_Matrix.dir/src/Matrix.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/src/Matrix.cpp -o CMakeFiles/calculate_A_Matrix.dir/src/Matrix.cpp.s

# Object files for target calculate_A_Matrix
calculate_A_Matrix_OBJECTS = \
"CMakeFiles/calculate_A_Matrix.dir/src/main.cpp.o" \
"CMakeFiles/calculate_A_Matrix.dir/src/helpers.cpp.o" \
"CMakeFiles/calculate_A_Matrix.dir/src/Matrix.cpp.o"

# External object files for target calculate_A_Matrix
calculate_A_Matrix_EXTERNAL_OBJECTS =

calculate_A_Matrix: CMakeFiles/calculate_A_Matrix.dir/src/main.cpp.o
calculate_A_Matrix: CMakeFiles/calculate_A_Matrix.dir/src/helpers.cpp.o
calculate_A_Matrix: CMakeFiles/calculate_A_Matrix.dir/src/Matrix.cpp.o
calculate_A_Matrix: CMakeFiles/calculate_A_Matrix.dir/build.make
calculate_A_Matrix: /usr/lib/x86_64-linux-gnu/libdolfin.so.2019.2.0.dev0
calculate_A_Matrix: /usr/lib/x86_64-linux-gnu/libpugixml.a
calculate_A_Matrix: /usr/lib/x86_64-linux-gnu/libboost_timer.so
calculate_A_Matrix: /usr/lib/x86_64-linux-gnu/libboost_chrono.so
calculate_A_Matrix: /usr/lib/x86_64-linux-gnu/hdf5/openmpi/libhdf5.so
calculate_A_Matrix: /usr/lib/x86_64-linux-gnu/libsz.so
calculate_A_Matrix: /usr/lib/x86_64-linux-gnu/libz.so
calculate_A_Matrix: /usr/lib/x86_64-linux-gnu/libdl.so
calculate_A_Matrix: /usr/lib/x86_64-linux-gnu/libm.so
calculate_A_Matrix: /usr/lib/slepcdir/slepc3.12/x86_64-linux-gnu-real/lib/libslepc_real.so
calculate_A_Matrix: /usr/lib/petscdir/petsc3.12/x86_64-linux-gnu-real/lib/libpetsc_real.so
calculate_A_Matrix: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi_cxx.so
calculate_A_Matrix: /usr/lib/x86_64-linux-gnu/openmpi/lib/libmpi.so
calculate_A_Matrix: CMakeFiles/calculate_A_Matrix.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable calculate_A_Matrix"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/calculate_A_Matrix.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/calculate_A_Matrix.dir/build: calculate_A_Matrix

.PHONY : CMakeFiles/calculate_A_Matrix.dir/build

CMakeFiles/calculate_A_Matrix.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/calculate_A_Matrix.dir/cmake_clean.cmake
.PHONY : CMakeFiles/calculate_A_Matrix.dir/clean

CMakeFiles/calculate_A_Matrix.dir/depend:
	cd /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/build /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/build /home/jiangl3/PODTherm-GP_Project_XML/Calculate_A/build/CMakeFiles/calculate_A_Matrix.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/calculate_A_Matrix.dir/depend

