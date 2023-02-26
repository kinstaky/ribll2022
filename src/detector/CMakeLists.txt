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
	PUBLIC ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
)


# tof library
add_library(tof tof.cpp ${DETECTOR_INCLUDE_DIR}/tof.h)
target_link_libraries(tof PUBLIC detector tof_event)


# dssd library
add_library(dssd dssd.cpp ${DETECTOR_INCLUDE_DIR}/dssd.h)
target_link_libraries(dssd PUBLIC detector dssd_event)


# # ppac library
# add_library(
# 	ppac
# 	ppac.cpp ${DETECTOR_INCLUDE_DIR}/ppac.h ${PROJECT_INCLUDE_DIR}/defs.h
# )
# target_include_directories(
# 	ppac
# 	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
# )
# target_link_libraries(
# 	ppac
# 	PUBLIC ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist
# )