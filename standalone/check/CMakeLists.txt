set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

# check chain
add_executable(check_chain check_chain.cpp)
target_link_libraries(check_chain PRIVATE chain_check)

# check pid
add_executable(check_t0_pid check_t0_pid.cpp)
target_link_libraries(check_t0_pid PRIVATE event t0)

# check three body information
add_executable(check_threebody_info check_threebody_info.cpp)
target_link_libraries(check_threebody_info PRIVATE event)

# check three body time
add_executable(check_threebody_time check_threebody_time.cpp)
target_link_libraries(check_threebody_time PRIVATE event ROOT::Graf)

# check angle calculation
add_executable(check_angle check_angle.cpp)
target_link_libraries(check_angle PRIVATE event)

# check TAFCsI shape
add_executable(check_tafcsi_shape check_tafcsi_shape.cpp)
target_include_directories(check_tafcsi_shape PRIVATE ROOT_INCLUDE_DIRS)
target_link_libraries(check_tafcsi_shape PRIVATE ROOT::Core ROOT::Graf)