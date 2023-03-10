set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

# check chain
add_executable(check_chain check_chain.cpp)
target_link_libraries(check_chain PRIVATE chain_check)