add_executable(taf_pid_image taf_pid_image.cpp)
target_link_libraries(
	taf_pid_image PRIVATE
	event ROOT::RIO ROOT::Tree ROOT::Graf ROOT::Gpad)

add_executable(taf_edge_q_thickness taf_edge_q_thickness.cpp)
target_link_libraries(
	taf_edge_q_thickness PRIVATE
	channel_v2_event csi_energy_calculator ROOT::Graf ROOT::Gpad
)