# show T0D2 shift
add_executable(show_d2_shift show_d2_shift.cpp)
target_link_libraries(show_d2_shift PRIVATE event ROOT::Core ROOT::RIO ROOT::Graf ROOT::Gpad)

# show T0D2 pixel
add_executable(show_d2_pixel show_d2_pixel.cpp)
target_link_libraries(show_d2_pixel PRIVATE event ROOT::Core ROOT::RIO ROOT::Graf ROOT::Gpad)

# correct hole flag
add_executable(correct_hole_flag correct_hole_flag.cpp)
target_include_directories(correct_hole_flag PRIVATE ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR)
target_link_libraries(correct_hole_flag PRIVATE ROOT::Core ROOT::RIO ROOT::Graf)

# hole model
add_executable(hole_model hole_model.cpp)
target_link_libraries(
	hole_model PRIVATE lost_energy_calculator ROOT::Core ROOT::RIO ROOT::Tree
)

# D1D2 energy relation
add_executable(hole_pid hole_pid.cpp)
target_link_libraries(
	hole_pid PRIVATE
	lost_energy_calculator d2_energy_calculator
	ROOT::Core ROOT::RIO ROOT::Tree
)

# pick events
add_executable(pick_hole_events pick_hole_events.cpp)
target_link_libraries(
	pick_hole_events PRIVATE event d2_energy_calculator ROOT::Core ROOT::RIO 
)
