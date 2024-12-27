add_executable(taf_pid_image taf_pid_image.cpp)
target_link_libraries(
	taf_pid_image PRIVATE
	event ROOT::RIO ROOT::Tree ROOT::Graf ROOT::Gpad)

add_executable(t0_pid_with_deutron t0_pid_with_deutron.cpp)
target_link_libraries(
	t0_pid_with_deutron PRIVATE
	event ROOT::RIO ROOT::Graf ROOT::Gpad
)

add_executable(taf_pid_with_t0 taf_pid_with_t0.cpp)
target_link_libraries(
	taf_pid_with_t0 PRIVATE
	event ROOT::RIO ROOT::Graf ROOT::Gpad
)

add_executable(taf_edge_q_thickness taf_edge_q_thickness.cpp)
target_link_libraries(
	taf_edge_q_thickness PRIVATE
	channel_v2_event csi_energy_calculator ROOT::Graf ROOT::Gpad
)

add_executable(spectrum_v3_q_range spectrum_v3_q_range.cpp)
target_include_directories(spectrum_v3_q_range PRIVATE ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(
	spectrum_v3_q_range PRIVATE
	ROOT::Graf ROOT::Gpad
)

add_executable(q_spectrum q_spectrum.cpp)
target_include_directories(q_spectrum PRIVATE ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(q_spectrum PRIVATE ROOT::Graf ROOT::Core ROOT::RIO ROOT::Tree ROOT::Gpad)