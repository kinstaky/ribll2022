set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

# check tabcsi mapping
add_executable(
	check_tacsi_mapping
	check_tacsi_mapping.cpp ${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	check_tacsi_mapping
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR
)
target_link_libraries(
	check_tacsi_mapping
	PUBLIC ROOT::Core ROOT::Tree ROOT::RIO ROOT::Hist ROOT::Gpad
)


# check chain
add_executable(check_chain check_chain.cpp)
target_link_libraries(check_chain PRIVATE chain_check)