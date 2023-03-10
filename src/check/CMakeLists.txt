# include directory
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(CHECK_INCLUDE_DIR ${PROJECT_INCLUDE_DIR}/check)

# check check library
add_library(chain_check chain_check.cpp ${CHECK_INCLUDE_DIR}/chain_check.h)
target_include_directories(chain_check PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(chain_check PUBLIC statistics)