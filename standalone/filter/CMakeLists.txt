# filter normalize
add_executable(filter_normalize filter_normalize.cpp)
target_link_libraries(filter_normalize PRIVATE event ROOT::Graf)

# t0d1d2 4He tail filter
add_executable(filter_he_tail filter_he_tail.cpp)
target_link_libraries(filter_he_tail PRIVATE event ROOT::Graf)

# t0d1d2 4He filter
add_executable(filter_he filter_he.cpp)
target_link_libraries(filter_he PRIVATE event)

# t0d1d2 line filter
add_executable(filter_t0d1d2_line filter_t0d1d2_line.cpp)
target_link_libraries(filter_t0d1d2_line PRIVATE event)

# time energy curve
add_executable(filter_time_curve filter_time_curve.cpp)
target_link_libraries(filter_time_curve PRIVATE detectors)

# filter dead range xia trigger
add_executable(filter_dead_range_xt filter_dead_range_xt.cpp)
target_link_libraries(filter_dead_range_xt PRIVATE event ROOT::Core)