# include directory
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(DETECTOR_INCLUDE_DIR ${PROJECT_INCLUDE_DIR}/detector)

# ppac library
add_library(
	ppac
	ppac.cpp ${DETECTOR_INCLUDE_DIR}/ppac.h ${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	ppac
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(
	ppac
	PUBLIC ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
)