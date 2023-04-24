# include directory
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(STATISTICS_INCLUDE_DIR ${PROJECT_INCLUDE_DIR}/statistics)

set(
	STATISTICS_HEADER_FILES
	${PROJECT_INCLUDE_DIR}/defs.h
	${STATISTICS_INCLUDE_DIR}/csv_line_reader.h
	${STATISTICS_INCLUDE_DIR}/statistics.h
	${STATISTICS_INCLUDE_DIR}/map_statistics.h
	${STATISTICS_INCLUDE_DIR}/align_statistics.h
	${STATISTICS_INCLUDE_DIR}/match_trigger_statistics.h
	${STATISTICS_INCLUDE_DIR}/normalize_statistics.h
	${STATISTICS_INCLUDE_DIR}/merge_case_statistics.h
	${STATISTICS_INCLUDE_DIR}/merge_statistics.h
	${STATISTICS_INCLUDE_DIR}/track_statistics.h
	${STATISTICS_INCLUDE_DIR}/track_dssd_statistics.h

	${STATISTICS_INCLUDE_DIR}/xt_period_statistics.h
	${STATISTICS_INCLUDE_DIR}/beam_identify_statistics.h
	${STATISTICS_INCLUDE_DIR}/t0_hit_statistics.h
)

set(
	STATISTICS_SOURCE_FILES
	csv_line_reader.cpp
	statistics.cpp
	map_statistics.cpp
	align_statistics.cpp
	match_trigger_statistics.cpp
	normalize_statistics.cpp
	merge_case_statistics.cpp
	merge_statistics.cpp
	track_statistics.cpp
	track_dssd_statistics.cpp

	xt_period_statistics.cpp
	beam_identify_statistics.cpp
	t0_hit_statistics.cpp
)

# statistics library
add_library(
	statistics
	${STATISTICS_SOURCE_FILES} ${STATISTICS_HEADER_FILES}
)
target_include_directories(statistics PUBLIC ${PROJECT_SOURCE_DIR})