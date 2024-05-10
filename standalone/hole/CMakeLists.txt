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
