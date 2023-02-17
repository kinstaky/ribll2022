# include directory
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

add_subdirectory(detector)
add_subdirectory(event)

# # defs library
# add_library(
# 	defs
# 	defs.cpp ${PROJECT_INCLUDE_DIR}/defs.h
# )
# target_include_directories(
# 	defs
# 	PUBLIC ${PROJECT_SOURCE_DIR}
# )

# mapper library
add_library(
	mapper
	mapper.cpp ${PROJECT_INCLUDE_DIR}/mapper.h ${PROJECT_INCLUDE_DIR}/defs.h
	${PROJECT_INCLUDE_DIR}/event/event.h
)
target_include_directories(
	mapper
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(
	mapper
	PUBLIC ROOT::Core ROOT::RIO ROOT::Tree
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
	PUBLIC ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
)

# telescope
add_library(
	telescope
	telescope.cpp ${PROJECT_INCLUDE_DIR}/telescope.h ${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	telescope
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(
	telescope
	PUBLIC detector ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
)