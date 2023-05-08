# filter normalize
add_executable(filter_normalize filter_normalize.cpp)
target_link_libraries(filter_normalize PRIVATE detectors)

# t0d1d2 line filter
add_executable(filter_t0d1d2_line filter_t0d1d2_line.cpp)
target_link_libraries(filter_t0d1d2_line PRIVATE event)