# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ferrari/work

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ferrari/work/build

# Include any dependencies generated for this target.
include CMakeFiles/pATTPClibSP.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/pATTPClibSP.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/pATTPClibSP.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/pATTPClibSP.dir/flags.make

CMakeFiles/pATTPClibSP.dir/src/mlesac.cc.o: CMakeFiles/pATTPClibSP.dir/flags.make
CMakeFiles/pATTPClibSP.dir/src/mlesac.cc.o: ../src/mlesac.cc
CMakeFiles/pATTPClibSP.dir/src/mlesac.cc.o: CMakeFiles/pATTPClibSP.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ferrari/work/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/pATTPClibSP.dir/src/mlesac.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/pATTPClibSP.dir/src/mlesac.cc.o -MF CMakeFiles/pATTPClibSP.dir/src/mlesac.cc.o.d -o CMakeFiles/pATTPClibSP.dir/src/mlesac.cc.o -c /home/ferrari/work/src/mlesac.cc

CMakeFiles/pATTPClibSP.dir/src/mlesac.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pATTPClibSP.dir/src/mlesac.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ferrari/work/src/mlesac.cc > CMakeFiles/pATTPClibSP.dir/src/mlesac.cc.i

CMakeFiles/pATTPClibSP.dir/src/mlesac.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pATTPClibSP.dir/src/mlesac.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ferrari/work/src/mlesac.cc -o CMakeFiles/pATTPClibSP.dir/src/mlesac.cc.s

CMakeFiles/pATTPClibSP.dir/src/ransac.cc.o: CMakeFiles/pATTPClibSP.dir/flags.make
CMakeFiles/pATTPClibSP.dir/src/ransac.cc.o: ../src/ransac.cc
CMakeFiles/pATTPClibSP.dir/src/ransac.cc.o: CMakeFiles/pATTPClibSP.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ferrari/work/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/pATTPClibSP.dir/src/ransac.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/pATTPClibSP.dir/src/ransac.cc.o -MF CMakeFiles/pATTPClibSP.dir/src/ransac.cc.o.d -o CMakeFiles/pATTPClibSP.dir/src/ransac.cc.o -c /home/ferrari/work/src/ransac.cc

CMakeFiles/pATTPClibSP.dir/src/ransac.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/pATTPClibSP.dir/src/ransac.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ferrari/work/src/ransac.cc > CMakeFiles/pATTPClibSP.dir/src/ransac.cc.i

CMakeFiles/pATTPClibSP.dir/src/ransac.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/pATTPClibSP.dir/src/ransac.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ferrari/work/src/ransac.cc -o CMakeFiles/pATTPClibSP.dir/src/ransac.cc.s

# Object files for target pATTPClibSP
pATTPClibSP_OBJECTS = \
"CMakeFiles/pATTPClibSP.dir/src/mlesac.cc.o" \
"CMakeFiles/pATTPClibSP.dir/src/ransac.cc.o"

# External object files for target pATTPClibSP
pATTPClibSP_EXTERNAL_OBJECTS =

libpATTPClibSP.so: CMakeFiles/pATTPClibSP.dir/src/mlesac.cc.o
libpATTPClibSP.so: CMakeFiles/pATTPClibSP.dir/src/ransac.cc.o
libpATTPClibSP.so: CMakeFiles/pATTPClibSP.dir/build.make
libpATTPClibSP.so: CMakeFiles/pATTPClibSP.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ferrari/work/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX shared library libpATTPClibSP.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/pATTPClibSP.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/pATTPClibSP.dir/build: libpATTPClibSP.so
.PHONY : CMakeFiles/pATTPClibSP.dir/build

CMakeFiles/pATTPClibSP.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/pATTPClibSP.dir/cmake_clean.cmake
.PHONY : CMakeFiles/pATTPClibSP.dir/clean

CMakeFiles/pATTPClibSP.dir/depend:
	cd /home/ferrari/work/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ferrari/work /home/ferrari/work /home/ferrari/work/build /home/ferrari/work/build /home/ferrari/work/build/CMakeFiles/pATTPClibSP.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/pATTPClibSP.dir/depend

