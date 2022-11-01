# include directory
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

# mapper library
add_library(
	mapper
	mapper.cpp ${PROJECT_INCLUDE_DIR}/mapper.h
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
	alignment.cpp ${PROJECT_INCLUDE_DIR}/alignment.h
)
target_include_directories(
	alignment
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(
	alignment
	PUBLIC ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
)