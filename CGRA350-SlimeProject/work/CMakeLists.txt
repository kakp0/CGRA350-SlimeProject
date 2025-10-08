################################################################
#                          Created by:                         #
#                    Ben Allen & Joshua Scott                  #
################################################################

cmake_minimum_required(VERSION 3.15)

# Include Extra CMake capabilities
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
include(CGRAFunctions)

# Project Name
set(CGRA_PROJECT "base" CACHE STRING "CGRA Project Name")
project("CGRA_PROJECT_${CGRA_PROJECT}" CXX C)

# Enable IDE Project Folders
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT base)


#########################################################
# Force Output Directories
#########################################################

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")



#########################################################
# Build Dependencies from Local 'ext' Folder
#########################################################
find_package(OpenGL REQUIRED)
find_package(OpenMP)

# Add the subdirectories containing the source code for your libraries.
# This tells CMake to go into these folders and build them.
add_subdirectory(ext/glfw)
add_subdirectory(ext/glew-1.10.0)
add_subdirectory(ext/glm)
add_subdirectory(ext/imgui)
add_subdirectory(ext/stb)


if (OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	add_definitions(-DCGRA_HAVE_OPENMP)
endif()



#########################################################
# Set Compiler Flags
#########################################################

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	add_compile_options(/std:c++latest)
	add_compile_options(/utf-8)
	add_compile_options(/W4 /MP)
	add_compile_options(/wd4800)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	add_compile_options("$<$<NOT:$<CONFIG:Debug>>:-O2>")
	add_compile_options(-std=c++1z -Wall -Wextra -pedantic)
	add_compile_options(-fvisibility=hidden)
	add_compile_options(-pthread -msse2)
	add_compile_options(-Werror=return-type)
	execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
	if (GCC_VERSION VERSION_GREATER 4.9 OR GCC_VERSION VERSION_EQUAL 4.9)
		add_compile_options(-fdiagnostics-color)
	endif()
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "^(Apple)?Clang$")
	add_compile_options("$<$<NOT:$<CONFIG:Debug>>:-O2>")
	add_compile_options(-std=c++1z -Wall -Wextra -pedantic)
	add_compile_options(-fvisibility=hidden)
	add_compile_options(-pthread -msse2)
	add_compile_options(-Werror=return-type)
endif()



#########################################################
# Source Files
#########################################################

add_subdirectory(src) 
add_subdirectory(res) 
set_property(TARGET ${CGRA_PROJECT} PROPERTY FOLDER "CGRA")

