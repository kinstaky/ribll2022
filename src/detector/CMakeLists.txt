# include directory
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(DETECTOR_INCLUDE_DIR ${PROJECT_INCLUDE_DIR}/detector)

# detector library
add_library(
	detector
	detector.cpp ${DETECTOR_INCLUDE_DIR}/detector.h ${PROJECT_INCLUDE_DIR}/defs.h
)
target_include_directories(
	detector
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(
	detector
	PUBLIC statistics ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
)


# tof library
add_library(tof tof.cpp ${DETECTOR_INCLUDE_DIR}/tof.h)
target_link_libraries(tof PUBLIC detector tof_event)


# vme trigger library
add_library(vme_trigger vme_trigger.cpp ${DETECTOR_INCLUDE_DIR}/vme_trigger.h)
target_link_libraries(
	vme_trigger
	PUBLIC detector trigger_event tof_event ppac_event dssd_event
)

# dssd library
add_library(dssd dssd.cpp ${DETECTOR_INCLUDE_DIR}/dssd.h)
target_link_libraries(dssd PUBLIC detector dssd_event)


# adssd library
add_library(adssd adssd.cpp ${DETECTOR_INCLUDE_DIR}/adssd.h)
target_link_libraries(adssd PUBLIC detector dssd_event)


#csi library
add_library(csi csi.cpp ${DETECTOR_INCLUDE_DIR}/csi.h)
target_link_libraries(csi PUBLIC detector csi_event)


# ppac library
add_library(ppac ppac.cpp ${DETECTOR_INCLUDE_DIR}/ppac.h)
target_link_libraries(ppac PUBLIC detector ppac_event)