# include directory
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(TELESCOPE_INCLUDE_DIR ${PROJECT_INCLUDE_DIR}/telescope)

# telescope library
add_library(
	telescope
	telescope.cpp ${TELESCOPE_INCLUDE_DIR}/telescope.h
	${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	telescope
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(
	telescope
	PUBLIC statistics ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
)

# taf library
add_library(
	taf
	taf.cpp ${TELESCOPE_INCLUDE_DIR}/telescope.h
)
target_link_libraries(taf PUBLIC telescope event)