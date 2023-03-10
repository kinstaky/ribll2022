# include directory
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

# single crate map
add_executable(smap single_map.cpp)
target_link_libraries(smap PUBLIC mapper)

# align
add_executable(align align.cpp)
target_link_libraries(align PUBLIC alignment)

# map and align
add_executable(map map.cpp)
target_link_libraries(map PUBLIC mapper alignment)

# match trigger
add_executable(match_trigger match_trigger.cpp)
target_link_libraries(match_trigger PUBLIC detectors)

# # correlation
# add_executable(correlate correlate.cpp)
# target_link_libraries(correlate PUBLIC detector ppac)

# # merge
# add_executable(merge merge.cpp)
# target_link_libraries(merge PUBLIC detector)

# # normalize
# add_executable(normalize normalize.cpp)
# target_link_libraries(normalize PUBLIC detector)

# # telescope
# add_executable(tele telescope.cpp)
# target_link_libraries(tele PUBLIC telescope)

# # ppac correlation
# add_executable(ppac_correlation ppac_correlation.cpp)
# target_include_directories(ppac_correlation PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR)
# target_link_libraries(ppac_correlation PUBLIC ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist)

# # ppac tracking
# add_executable(track track.cpp)
# target_link_libraries(track PUBLIC ppac)

# # reftof
# add_executable(reftof reftof.cpp)
# target_include_directories(reftof PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR)
# target_link_libraries(reftof PUBLIC ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist)

add_subdirectory(check)
add_subdirectory(show)