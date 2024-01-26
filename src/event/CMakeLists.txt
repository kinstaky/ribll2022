# include directory
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(EVENT_INCLUDE_DIR ${PROJECT_INCLUDE_DIR}/event)

set(
	EVENT_HEADER_FILES
	${PROJECT_INCLUDE_DIR}/defs.h
	${EVENT_INCLUDE_DIR}/event.h
	${EVENT_INCLUDE_DIR}/csi_event.h
	${EVENT_INCLUDE_DIR}/dssd_event.h
	${EVENT_INCLUDE_DIR}/ppac_event.h
	${EVENT_INCLUDE_DIR}/ssd_event.h
	${EVENT_INCLUDE_DIR}/tof_event.h
	${EVENT_INCLUDE_DIR}/trigger_event.h
	${EVENT_INCLUDE_DIR}/ta_event.h
	${EVENT_INCLUDE_DIR}/t0_event.h
	${EVENT_INCLUDE_DIR}/particle_type_event.h
	${EVENT_INCLUDE_DIR}/particle_event.h
	${EVENT_INCLUDE_DIR}/channel_event.h
	${EVENT_INCLUDE_DIR}/generate_event.h
	${EVENT_INCLUDE_DIR}/detect_event.h
	${EVENT_INCLUDE_DIR}/filter_event.h
	${EVENT_INCLUDE_DIR}/pd_info_event.h
	${EVENT_INCLUDE_DIR}/threebody_info_event.h
)

set(
	EVENT_SOURCE_FILES
	csi_event.cpp
	dssd_event.cpp
	ppac_event.cpp
	ssd_event.cpp
	tof_event.cpp
	trigger_event.cpp
	ta_event.cpp
	t0_event.cpp
	particle_type_event.cpp
	particle_event.cpp
	channel_event.cpp
	generate_event.cpp
	detect_event.cpp
	filter_event.cpp
	pd_info_event.cpp
	threebody_info_event.cpp
)

# event library
add_library(
	event
	${EVENT_SOURCE_FILES} ${EVENT_HEADER_FILES}
)
target_include_directories(
	event
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR
)
target_link_libraries(event PUBLIC ROOT::Tree)