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
CMAKE_SOURCE_DIR = /mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src

# Include any dependencies generated for this target.
include CMakeFiles/oct2grid.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/oct2grid.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/oct2grid.dir/flags.make

CMakeFiles/oct2grid.dir/oct2grid.c.o: CMakeFiles/oct2grid.dir/flags.make
CMakeFiles/oct2grid.dir/oct2grid.c.o: oct2grid.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/oct2grid.dir/oct2grid.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/oct2grid.dir/oct2grid.c.o   -c /mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/oct2grid.c

CMakeFiles/oct2grid.dir/oct2grid.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/oct2grid.dir/oct2grid.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/oct2grid.c > CMakeFiles/oct2grid.dir/oct2grid.c.i

CMakeFiles/oct2grid.dir/oct2grid.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/oct2grid.dir/oct2grid.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/oct2grid.c -o CMakeFiles/oct2grid.dir/oct2grid.c.s

# Object files for target oct2grid
oct2grid_OBJECTS = \
"CMakeFiles/oct2grid.dir/oct2grid.c.o"

# External object files for target oct2grid
oct2grid_EXTERNAL_OBJECTS = \
"/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles/GRID_LIB_OBJS.dir/GridLib.c.o" \
"/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles/GRID_LIB_OBJS.dir/util.c.o" \
"/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles/GRID_LIB_OBJS.dir/geo.c.o" \
"/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles/GRID_LIB_OBJS.dir/octtree/octtree.c.o" \
"/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles/GRID_LIB_OBJS.dir/io/json_io.c.o" \
"/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles/GRID_LIB_OBJS.dir/io/jReadWrite/source/jRead.c.o" \
"/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles/GRID_LIB_OBJS.dir/io/jReadWrite/source/jWrite.c.o" \
"/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles/GRID_LIB_OBJS.dir/alomax_matrix/alomax_matrix.c.o" \
"/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles/GRID_LIB_OBJS.dir/alomax_matrix/eigv.c.o" \
"/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles/GRID_LIB_OBJS.dir/alomax_matrix/alomax_matrix_svd.c.o" \
"/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles/GRID_LIB_OBJS.dir/matrix_statistics/matrix_statistics.c.o" \
"/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles/GRID_LIB_OBJS.dir/vector/vector.c.o" \
"/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles/GRID_LIB_OBJS.dir/ran1/ran1.c.o" \
"/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles/GRID_LIB_OBJS.dir/map_project.c.o"

bin/oct2grid: CMakeFiles/oct2grid.dir/oct2grid.c.o
bin/oct2grid: CMakeFiles/GRID_LIB_OBJS.dir/GridLib.c.o
bin/oct2grid: CMakeFiles/GRID_LIB_OBJS.dir/util.c.o
bin/oct2grid: CMakeFiles/GRID_LIB_OBJS.dir/geo.c.o
bin/oct2grid: CMakeFiles/GRID_LIB_OBJS.dir/octtree/octtree.c.o
bin/oct2grid: CMakeFiles/GRID_LIB_OBJS.dir/io/json_io.c.o
bin/oct2grid: CMakeFiles/GRID_LIB_OBJS.dir/io/jReadWrite/source/jRead.c.o
bin/oct2grid: CMakeFiles/GRID_LIB_OBJS.dir/io/jReadWrite/source/jWrite.c.o
bin/oct2grid: CMakeFiles/GRID_LIB_OBJS.dir/alomax_matrix/alomax_matrix.c.o
bin/oct2grid: CMakeFiles/GRID_LIB_OBJS.dir/alomax_matrix/eigv.c.o
bin/oct2grid: CMakeFiles/GRID_LIB_OBJS.dir/alomax_matrix/alomax_matrix_svd.c.o
bin/oct2grid: CMakeFiles/GRID_LIB_OBJS.dir/matrix_statistics/matrix_statistics.c.o
bin/oct2grid: CMakeFiles/GRID_LIB_OBJS.dir/vector/vector.c.o
bin/oct2grid: CMakeFiles/GRID_LIB_OBJS.dir/ran1/ran1.c.o
bin/oct2grid: CMakeFiles/GRID_LIB_OBJS.dir/map_project.c.o
bin/oct2grid: CMakeFiles/oct2grid.dir/build.make
bin/oct2grid: CMakeFiles/oct2grid.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable bin/oct2grid"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/oct2grid.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/oct2grid.dir/build: bin/oct2grid

.PHONY : CMakeFiles/oct2grid.dir/build

CMakeFiles/oct2grid.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/oct2grid.dir/cmake_clean.cmake
.PHONY : CMakeFiles/oct2grid.dir/clean

CMakeFiles/oct2grid.dir/depend:
	cd /mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src /mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src /mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src /mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src /mnt/e/Documents/Lavori_Tirocini/Assegno_CRS_22_23/TomoNE/python/modules/NllocPy/NllocPy/src/CMakeFiles/oct2grid.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/oct2grid.dir/depend

