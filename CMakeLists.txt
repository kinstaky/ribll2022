# Require the cmake minimum version.
cmake_minimum_required(VERSION 3.10)


# Set project properties.
project(
	Ribll2022
	VERSION 0.1.0
	DESCRIPTION "A project to analysis data from ribll 1 experiment in 2022 by PKU"
	LANGUAGES CXX
)


# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if (EXISTS "${LOC_PATH}")
	message(FATAL_ERROR "You cannot build in a source directory (or any directory with"
						"CmakeLists.txt file). Please make a build subdirectory. Feel free"
						"to remove CMakeCache.txt and CMakeFiles.")
endif()

# Set the build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()


# Only do these if this is the main project, and not if it is included
# through add_subdirectory
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

	# Let's ensure -std=c++xx instead of -std=g++xx
	set(CMAKE_CXX_EXTENSIONS OFF)

	# Set the default cxx standard
	set(CMAKE_CXX_STANDARD 17)

	# Set default compile options
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra")

	# Let's nicely support folders in IDE's
	set_property(GLOBAL PROPERTY USE_FOLDERS ON)

	# # Enable testing if it's in the main app.
	# include (CTest)


	# # Docs only available if this is the main app
	# find_package(Doxygen)
	# if(Doxygen_FOUND)
	# 	add_subdirectory(docs)
	# else()
	# 	message(STATUS "Doxygen not found, not building docs")
	# endif()

endif()

# find root
find_package(ROOT 6.26 REQUIRED COMPONENTS Core RIO Tree Hist Graf)

option(BUILD_ALL "Builds all parts" ON)
option(BUILD_ANALYSIS "Builds analysis part" OFF)
option(BUILD_NUMBERIC_OPTIMIZE "Builds numberic optmization part" OFF)
option(BUILD_SIMULATION "Builds simulation part" OFF)

if (BUILD_ALL)
	set(BUILD_ANALYSIS ON)
	set(BUILD_NUMBERIC_OPTIMIZE ON)
	set(BUILD_SIMULATION ON)
endif()

function(show_config_options)
	message("")
	message("Build Properties")
	message("  c++ compiler.................: " ${CMAKE_CXX_COMPILER})
	message("  c++ stardard.................: " ${CMAKE_CXX_STANDARD})
	message("  cmake........................: " ${CMAKE_VERSION})
	message("  cmake build type.............: " ${CMAKE_BUILD_TYPE})
	message("")
	message("Build Options")
	message("  BUILD_ALL....................: " ${BUILD_ALL})
	message("  BUILD_ANALYSIS...............: " ${BUILD_ANALYSIS})
	message("  BUILD_NUMBERIC_OPTIMIZE......: " ${BUILD_NUMBERIC_OPTIMIZE})
	message("  BUILD_SIMULATION.............: " ${BUILD_SIMULATION})
endfunction()

show_config_options()

if (BUILD_ANALYSIS)
	# The complied library code
	add_subdirectory(src)
	add_subdirectory(src/event)
	add_subdirectory(src/calculator)
	# Standalone executable program.
	add_subdirectory(standalone)
endif()

if (BUILD_NUMBERIC_OPTIMIZE)
	add_subdirectory(standalone/optimize)
	if (NOT BUILD_ANALYSIS)
		add_subdirectory(src/event)
		add_subdirectory(src/calculator)
	endif()
endif()

if (BUILD_SIMULATION)
	add_subdirectory(standalone/simulate)
	if (NOT BUILD_ANALYSIS AND NOT BUILD_NUMBERIC_OPTIMIZE)
		add_subdirectory(src/event)
		add_subdirectory(src/calculator)
	endif()
endif()