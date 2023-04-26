set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

# show tabcsi mapping
add_executable(
	show_tacsi_mapping
	show_tacsi_mapping.cpp ${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	show_tacsi_mapping
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR
)
target_link_libraries(
	show_tacsi_mapping
	PUBLIC ROOT::Core ROOT::Tree ROOT::RIO ROOT::Hist ROOT::Gpad
)

# show xt period
add_executable(
	show_xt_period
	show_xt_period.cpp ${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	show_xt_period
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR
)
target_link_libraries(
	show_xt_period
	PUBLIC statistics ROOT::Core ROOT::Tree ROOT::RIO ROOT::Hist
)

# show T0 hit
add_executable(show_t0_hit show_t0_hit.cpp ${PROJECT_INCLUDE_DIR}/defs.h)
target_include_directories(
	show_t0_hit
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR
)
target_link_libraries(
	show_t0_hit
	PUBLIC statistics ROOT::Core ROOT::Tree ROOT::RIO
)

# show gdc offset
add_executable(show_gdc_offset show_gdc_offset.cpp)
target_include_directories(
	show_gdc_offset
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR
)
target_link_libraries(
	show_gdc_offset
	PUBLIC ROOT::Core ROOT::Tree ROOT::RIO ROOT::Graf
)

# show align group
add_executable(show_align_group show_align_group.cpp)
target_include_directories(
	show_align_group
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR
)
target_link_libraries(
	show_align_group
	PUBLIC alignment
)

# show normalize result
add_executable(show_normalize show_normalize.cpp)
target_link_libraries(
	show_normalize PRIVATE detectors
)

# show t0 center
add_executable(show_t0_center show_t0_center.cpp)
target_link_libraries(
	show_t0_center PRIVATE event ROOT::Graf
)

# show T0 pid
add_executable(show_t0_pid show_t0_pid.cpp)
target_link_libraries(show_t0_pid PRIVATE event ROOT::Hist)

# show TAF pid
add_executable(show_taf_pid show_taf_pid.cpp)
target_link_libraries(show_taf_pid PRIVATE event ROOT::Hist)


# show PPAC center
add_executable(show_ppac_center show_ppac_center.cpp)
target_link_libraries(
	show_ppac_center PRIVATE event ROOT::Graf
)