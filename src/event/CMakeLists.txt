# include directory
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(EVENT_INCLUDE_DIR ${PROJECT_INCLUDE_DIR}/event)

# tof event library
add_library(
	tof_event
	tof_event.cpp ${EVENT_INCLUDE_DIR}/tof_event.h
	${EVENT_INCLUDE_DIR}/event.h ${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	tof_event
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(tof_event PUBLIC ROOT::Tree)

# dssd event library
add_library(
	dssd_event
	dssd_event.cpp ${EVENT_INCLUDE_DIR}/dssd_event.h
	${EVENT_INCLUDE_DIR}/event.h ${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	dssd_event
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(dssd_event PUBLIC ROOT::Tree)

# csi event library
add_library(
	csi_event
	csi_event.cpp ${EVENT_INCLUDE_DIR}/csi_event.h
	${EVENT_INCLUDE_DIR}/event.h ${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	csi_event
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(csi_event PUBLIC ROOT::Tree)

# ppac event library
add_library(
	ppac_event
	ppac_event.cpp ${EVENT_INCLUDE_DIR}/ppac_event.h
	${EVENT_INCLUDE_DIR}/event.h ${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	ppac_event
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIR
)
target_link_libraries(ppac_event PUBLIC ROOT::Tree)