# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

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
CMAKE_COMMAND = "/Users/peterzheng/Library/Application Support/JetBrains/Toolbox/apps/CLion/ch-0/172.4343.16/CLion.app/Contents/bin/cmake/bin/cmake"

# The command to remove a file.
RM = "/Users/peterzheng/Library/Application Support/JetBrains/Toolbox/apps/CLion/ch-0/172.4343.16/CLion.app/Contents/bin/cmake/bin/cmake" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/peterzheng/GitHubProjects/OOP_Matrix_2017Fall

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/peterzheng/GitHubProjects/OOP_Matrix_2017Fall/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/OOP_Matrix_2017Fall.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/OOP_Matrix_2017Fall.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/OOP_Matrix_2017Fall.dir/flags.make

CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.o: CMakeFiles/OOP_Matrix_2017Fall.dir/flags.make
CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.o: ../BasicTest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/peterzheng/GitHubProjects/OOP_Matrix_2017Fall/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.o"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.o -c /Users/peterzheng/GitHubProjects/OOP_Matrix_2017Fall/BasicTest.cpp

CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.i"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/peterzheng/GitHubProjects/OOP_Matrix_2017Fall/BasicTest.cpp > CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.i

CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.s"
	/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/peterzheng/GitHubProjects/OOP_Matrix_2017Fall/BasicTest.cpp -o CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.s

CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.o.requires:

.PHONY : CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.o.requires

CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.o.provides: CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.o.requires
	$(MAKE) -f CMakeFiles/OOP_Matrix_2017Fall.dir/build.make CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.o.provides.build
.PHONY : CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.o.provides

CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.o.provides.build: CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.o


# Object files for target OOP_Matrix_2017Fall
OOP_Matrix_2017Fall_OBJECTS = \
"CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.o"

# External object files for target OOP_Matrix_2017Fall
OOP_Matrix_2017Fall_EXTERNAL_OBJECTS =

OOP_Matrix_2017Fall: CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.o
OOP_Matrix_2017Fall: CMakeFiles/OOP_Matrix_2017Fall.dir/build.make
OOP_Matrix_2017Fall: CMakeFiles/OOP_Matrix_2017Fall.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/peterzheng/GitHubProjects/OOP_Matrix_2017Fall/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable OOP_Matrix_2017Fall"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/OOP_Matrix_2017Fall.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/OOP_Matrix_2017Fall.dir/build: OOP_Matrix_2017Fall

.PHONY : CMakeFiles/OOP_Matrix_2017Fall.dir/build

CMakeFiles/OOP_Matrix_2017Fall.dir/requires: CMakeFiles/OOP_Matrix_2017Fall.dir/BasicTest.cpp.o.requires

.PHONY : CMakeFiles/OOP_Matrix_2017Fall.dir/requires

CMakeFiles/OOP_Matrix_2017Fall.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/OOP_Matrix_2017Fall.dir/cmake_clean.cmake
.PHONY : CMakeFiles/OOP_Matrix_2017Fall.dir/clean

CMakeFiles/OOP_Matrix_2017Fall.dir/depend:
	cd /Users/peterzheng/GitHubProjects/OOP_Matrix_2017Fall/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/peterzheng/GitHubProjects/OOP_Matrix_2017Fall /Users/peterzheng/GitHubProjects/OOP_Matrix_2017Fall /Users/peterzheng/GitHubProjects/OOP_Matrix_2017Fall/cmake-build-debug /Users/peterzheng/GitHubProjects/OOP_Matrix_2017Fall/cmake-build-debug /Users/peterzheng/GitHubProjects/OOP_Matrix_2017Fall/cmake-build-debug/CMakeFiles/OOP_Matrix_2017Fall.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/OOP_Matrix_2017Fall.dir/depend

