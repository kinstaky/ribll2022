# simulate transfer reaction
add_executable(sim_transfer simulate_transfer_reaction.cpp)
target_include_directories(sim_transfer PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(sim_transfer PRIVATE ROOT::Core ROOT::Graf ROOT::RIO)

# generate simulate data
add_executable(generate generate.cpp)
target_link_libraries(generate PRIVATE event target_energy_calculator ROOT::Core ROOT::RIO ROOT::Graf)

# simulate detecting generated data
add_executable(detect detect.cpp)
target_link_libraries(detect PRIVATE event range_energy_calculator ROOT::Core ROOT::RIO ROOT::Graf)

# simulate rebuilding detected data
add_executable(sim_rebuild simulate_rebuild.cpp)
target_link_libraries(sim_rebuild PRIVATE event ROOT::RIO ROOT::Graf)

# track simulate data
add_executable(sim_track simulate_track.cpp)
target_link_libraries(sim_track PRIVATE event ROOT::Graf)

# single ppac track
add_executable(single_ppac_track single_ppac_track.cpp)
target_link_libraries(single_ppac_track PRIVATE event ROOT::Graf)

# resist strip
add_executable(sim_resist_strip simulate_resist_strip.cpp)
target_link_libraries(sim_resist_strip PRIVATE event ROOT::Graf)

# simulate live time
add_executable(live_time live_time.cpp)
target_include_directories(live_time PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(live_time PRIVATE ROOT::Core ROOT::RIO ROOT::Graf)