# include directory
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(MAPPER_INCLUDE_DIR ${PROJECT_INCLUDE_DIR}/mapper)

set(
	MAPPER_HEADER_FILES
	${PROJECT_INCLUDE_DIR}/defs.h
	${MAPPER_INCLUDE_DIR}/xia_mapper.h
	${MAPPER_INCLUDE_DIR}/crate0_mapper.h
	${MAPPER_INCLUDE_DIR}/crate1_mapper.h
	${MAPPER_INCLUDE_DIR}/crate2_mapper.h
	${MAPPER_INCLUDE_DIR}/crate3_mapper.h
	${MAPPER_INCLUDE_DIR}/trace_mapper.h
)

set(
	MAPPER_SOURCE_FILES
	xia_mapper.cpp
	crate0_mapper.cpp
	crate1_mapper.cpp
	crate2_mapper.cpp
	crate3_mapper.cpp
	trace_mapper.cpp
)

# MAPPER library
add_library(
	mapper
	${MAPPER_SOURCE_FILES} ${MAPPER_HEADER_FILES}
)
target_include_directories(
	mapper
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(
	mapper
	PUBLIC statistics event ROOT::Core ROOT::RIO ROOT::Tree
)