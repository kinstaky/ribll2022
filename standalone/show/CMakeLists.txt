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