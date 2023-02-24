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
find_package(ROOT 6.26 REQUIRED COMPONENTS Core RIO Tree Hist)

# The complied library code
add_subdirectory(src)

# Standalone executable program.
add_subdirectory(standalone)

# # add test
# if (BUILD_TESTING)
# 	set(INSTALL_GTEST OFF)
# 	add_subdirectory(test)
# endif()