# include directory
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

# mapping
add_executable(mapping mapping.cpp)
target_link_libraries(mapping PUBLIC mapper)

# align
add_executable(align align.cpp)
target_link_libraries(align PUBLIC alignment)