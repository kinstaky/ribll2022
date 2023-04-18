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
	PUBLIC
	statistics
	range_energy_calculator delta_energy_calculator csi_energy_calculator
	ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Graf
)

# taf library
add_library(
	taf
	taf.cpp ${TELESCOPE_INCLUDE_DIR}/taf.h
)
target_link_libraries(taf PUBLIC telescope event)

# t0 library
add_library(
	t0
	t0.cpp ${TELESCOPE_INCLUDE_DIR}/t0.h
)
target_link_libraries(t0 PUBLIC telescope event)