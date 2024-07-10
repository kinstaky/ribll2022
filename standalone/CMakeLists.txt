# include directory
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

# single crate map
add_executable(smap single_map.cpp)
target_link_libraries(smap PRIVATE mapper)

# align
add_executable(align align.cpp)
target_link_libraries(align PRIVATE alignment)

# map and align
add_executable(map map.cpp)
target_link_libraries(map PRIVATE mapper alignment)

# map trace
add_executable(tmap trace_map.cpp)
target_link_libraries(tmap PRIVATE mapper)

# match trigger
add_executable(match_trigger match_trigger.cpp)
target_link_libraries(match_trigger PRIVATE detectors)

# reference time
add_executable(ref_time ref_time.cpp)
target_link_libraries(ref_time time_reference)

# normalize time
add_executable(normalize_time normalize_time.cpp)
target_link_libraries(normalize_time PRIVATE detectors)

# fit time-energy curve
add_executable(fit_time_curve fit_time_curve.cpp)
target_link_libraries(fit_time_curve PRIVATE detectors)

# analyze time
add_executable(analyze_time analyze_time.cpp)
target_link_libraries(analyze_time PRIVATE detectors)

# beam identify
add_executable(beam_identify beam_identify.cpp)
target_link_libraries(beam_identify PRIVATE tof)

# normalize
add_executable(normalize normalize.cpp)
target_link_libraries(normalize PRIVATE detectors)

# normalize result
add_executable(normalize_result normalize_result.cpp)
target_link_libraries(normalize_result PRIVATE detectors)

# cut beam threshold
add_executable(cut_beam_threshold cut_beam_threshold.cpp)
target_link_libraries(cut_beam_threshold PRIVATE detectors)

# merge
add_executable(merge merge.cpp)
target_link_libraries(merge PRIVATE detectors)

# track
add_executable(track track.cpp)
target_link_libraries(track PRIVATE ppac telescopes)

# merge and track
add_executable(merge_track merge_and_track.cpp)
target_link_libraries(merge_track PRIVATE t0)

# merge binding events in T0D1 and T0D2
add_executable(merge_bind merge_bind.cpp)
target_link_libraries(merge_bind PRIVATE event dssd)

# analyze trace
add_executable(analyze_trace analyze_trace.cpp)
target_link_libraries(analyze_trace PRIVATE detectors telescopes)

# pid
add_executable(pid particle_identify.cpp)
target_link_libraries(pid PRIVATE telescopes)

# calibrate
add_executable(calibrate calibrate.cpp)
target_link_libraries(calibrate PRIVATE detectors telescopes)

# calibrate result
add_executable(calibrate_result calibrate_result.cpp)
target_link_libraries(calibrate_result PRIVATE telescopes)

# particle
add_executable(rebuild rebuild.cpp)
target_link_libraries(rebuild PRIVATE telescopes)

# coincide
add_executable(coincide coincide.cpp)
target_link_libraries(coincide PRIVATE channel)

# info-pd
add_executable(info_pd info_pd.cpp)
target_link_libraries(info_pd PRIVATE event ROOT::Graf)

# info-3body
add_executable(info_threebody info_threebody.cpp)
target_link_libraries(info_threebody PRIVATE event csi_energy_calculator d2_energy_calculator)

# # info-3body-extra
# add_executable(info_threebody_extra info_threebody_extra.cpp)
# target_link_libraries(info_threebody_extra PRIVATE event)

# compare 3body information between 9Be and 10Be
add_executable(compare_be_info compare_be_info.cpp)
target_link_libraries(compare_be_info PRIVATE event d2_energy_calculator)


add_subdirectory(calculate)
add_subdirectory(check)
add_subdirectory(filter)
add_subdirectory(show)
add_subdirectory(simulate)
add_subdirectory(summary)
add_subdirectory(optimize)
add_subdirectory(spectrum)
add_subdirectory(hole)