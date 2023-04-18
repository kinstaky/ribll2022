# include directory
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

add_subdirectory(check)
add_subdirectory(detector)
add_subdirectory(event)
add_subdirectory(mapper)
add_subdirectory(statistics)
add_subdirectory(telescope)
add_subdirectory(calculator)

# alignment
add_library(
	alignment
	alignment.cpp ${PROJECT_INCLUDE_DIR}/alignment.h ${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	alignment
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(
	alignment
	PUBLIC statistics ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
)

# detectors
add_library(detectors detectors.cpp ${PROJECT_INCLUDE_DIR}/detectors.h)
target_link_libraries(detectors PUBLIC adssd csi dssd ppac ssd tof vme_trigger)

# time reference
add_library(
	time_reference
	time_reference.cpp ${PROJECT_INCLUDE_DIR}/time_reference.h
)
target_include_directories(
	time_reference
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(
	time_reference
	PUBLIC event ROOT::Core ROOT::RIO ROOT::Tree
)

# telescopes
add_library(
	telescopes
	telescopes.cpp ${PROJECT_INCLUDE_DIR}/telescopes.h
)
target_link_libraries(
	telescopes
	PUBLIC taf t0
)