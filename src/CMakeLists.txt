# include directory
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

add_subdirectory(check)
add_subdirectory(detector)
add_subdirectory(event)

# statistics library
add_library(
	statistics
	statistics.cpp
	${PROJECT_INCLUDE_DIR}/statistics.h
	${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	statistics PUBLIC ${PROJECT_SOURCE_DIR}
)

# mapper library
add_library(
	mapper
	mapper.cpp ${PROJECT_INCLUDE_DIR}/mapper.h ${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	mapper
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(
	mapper
	PUBLIC statistics ppac_event dssd_event ROOT::Core ROOT::RIO ROOT::Tree
)

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
target_link_libraries(detectors PUBLIC tof dssd adssd csi vme_trigger)

# # telescope
# add_library(
# 	telescope
# 	telescope.cpp ${PROJECT_INCLUDE_DIR}/telescope.h ${PROJECT_INCLUDE_DIR}/defs.h
# )
# target_include_directories(
# 	telescope
# 	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
# )
# target_link_libraries(
# 	telescope
# 	PUBLIC detector ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
# )