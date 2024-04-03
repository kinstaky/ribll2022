set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)

# show tabcsi mapping
add_executable(
	show_tacsi_mapping
	show_tacsi_mapping.cpp ${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	show_tacsi_mapping
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR
)
target_link_libraries(
	show_tacsi_mapping
	PUBLIC ROOT::Core ROOT::Tree ROOT::RIO ROOT::Hist ROOT::Gpad
)

# show xt period
add_executable(
	show_event_interval
	show_event_interval.cpp ${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	show_event_interval
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR
)
target_link_libraries(
	show_event_interval
	PUBLIC statistics ROOT::Core ROOT::Tree ROOT::RIO ROOT::Hist
)

# show T0 hit
add_executable(show_t0_hit show_t0_hit.cpp ${PROJECT_INCLUDE_DIR}/defs.h)
target_include_directories(
	show_t0_hit
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR
)
target_link_libraries(
	show_t0_hit
	PUBLIC statistics ROOT::Core ROOT::Tree ROOT::RIO
)

# show gdc offset
add_executable(show_gdc_offset show_gdc_offset.cpp)
target_include_directories(
	show_gdc_offset
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR
)
target_link_libraries(
	show_gdc_offset
	PUBLIC ROOT::Core ROOT::Tree ROOT::RIO ROOT::Graf
)

# show align group
add_executable(show_align_group show_align_group.cpp)
target_include_directories(
	show_align_group
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR
)
target_link_libraries(
	show_align_group
	PUBLIC alignment
)

# show normalize result
add_executable(show_normalize show_normalize.cpp)
target_link_libraries(
	show_normalize PRIVATE detectors
)

# show T0 center
add_executable(show_t0_center show_t0_center.cpp)
target_link_libraries(
	show_t0_center PRIVATE event statistics ROOT::Graf
)

# show T0 inverse tracking
add_executable(show_t0_inverse_tracking show_t0_inverse_tracking.cpp)
target_link_libraries(
	show_t0_inverse_tracking PRIVATE event ROOT::Graf
)

# show T0 pid
add_executable(show_t0_pid show_t0_pid.cpp)
target_link_libraries(show_t0_pid PRIVATE event ROOT::Hist)

# show TAF pid
add_executable(show_taf_pid show_taf_pid.cpp)
target_link_libraries(show_taf_pid PRIVATE event ROOT::Hist)

# show PPAC center
add_executable(show_ppac_center show_ppac_center.cpp)
target_link_libraries(
	show_ppac_center PRIVATE event statistics ROOT::Graf
)

# show PPAC sum
add_executable(show_ppac_change show_ppac_change.cpp)
target_include_directories(
	show_ppac_change PRIVATE ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR
)
target_link_libraries(
	show_ppac_change PRIVATE ROOT::Core ROOT::RIO ROOT::Graf
)

# show PPAC effect
add_executable(show_ppac_effect show_ppac_effect.cpp)
target_link_libraries(
	show_ppac_effect PRIVATE event ROOT::Core ROOT::RIO ROOT::Graf
)

# show pd reaction single PPAC track
add_executable(show_pd_single_ppac show_pd_single_ppac.cpp)
target_link_libraries(
	show_pd_single_ppac PRIVATE event ROOT::Core ROOT::RIO ROOT::Graf
)

# show time-energy
add_executable(show_time_energy show_time_energy.cpp)
target_link_libraries(show_time_energy PRIVATE detectors event ROOT::Hist)

# show time-energy hit
add_executable(show_time_energy_hit show_time_energy_hit.cpp)
target_link_libraries(show_time_energy_hit PRIVATE detectors event ROOT::Hist)

# show t0 calibration result
add_executable(show_calibration show_calibration.cpp)
target_link_libraries(show_calibration PRIVATE t0)

# show pid width
add_executable(show_pid_width show_pid_width.cpp)
target_link_libraries(show_pid_width PRIVATE t0)

# show T0D3 hit in tail without tracking
add_executable(show_d3_hit show_d3_hit.cpp)
target_link_libraries(show_d3_hit PRIVATE event ROOT::Core ROOT::Graf)

# show T0D2 shift
add_executable(show_d2_shift show_d2_shift.cpp)
target_link_libraries(show_d2_shift PRIVATE event ROOT::Core ROOT::RIO ROOT::Graf ROOT::Gpad)

# show T0D2 pixel
add_executable(show_d2_pixel show_d2_pixel.cpp)
target_link_libraries(show_d2_pixel PRIVATE event ROOT::Core ROOT::RIO ROOT::Graf ROOT::Gpad)

# show resist strip
add_executable(show_resist_strip show_resist_strip.cpp)
target_link_libraries(show_resist_strip PRIVATE event ROOT::Graf)

# show tafd deduce thickness
add_executable(show_tafd_deduce_thick show_tafd_deduce_thickness.cpp)
target_link_libraries(
	show_tafd_deduce_thick PRIVATE
	event range_energy_calculator ROOT::Core ROOT::RIO ROOT::Graf
)

# show merge counts
add_executable(show_merge_counts show_merge_counts.cpp)
target_link_libraries(
	show_merge_counts PRIVATE
	event ROOT::Core ROOT::RIO ROOT::Graf
)

# show t0 rebuild rate
add_executable(show_t0_rebuild_rate show_t0_rebuild_rate.cpp)
target_link_libraries(
	show_t0_rebuild_rate PRIVATE
	event ROOT::Core ROOT::Graf
)

# show calculate T0D2 energy
add_executable(show_calculate_t0d2 show_calculate_t0d2.cpp)
target_link_libraries(
	show_calculate_t0d2 PRIVATE event d2_energy_calculator
)

# show T0 PID with angle
add_executable(show_t0_pid_angle show_t0_pid_angle.cpp)
target_link_libraries(
	show_t0_pid_angle PRIVATE event ROOT::Core ROOT::RIO ROOT::Graf
)