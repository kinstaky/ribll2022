# set(glog_DIR /data/d1/pwl/glog/lib64/cmake/glog/)
find_package(Ceres REQUIRED PATHS /data/d1/pwl/ceres-solver/)

# include directory
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

add_executable(threebody_opt threebody_optimize.cpp)
target_include_directories(
	threebody_opt PRIVATE ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(
	threebody_opt PRIVATE
	event delta_energy_calculator Ceres::ceres
)

add_executable(pd_group pd_group.cpp)
target_include_directories(
	pd_group PRIVATE ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(
	pd_group PRIVATE Ceres::ceres ROOT::Core ROOT::RIO ROOT::Tree ROOT::Graf
)

add_executable(threebody_group threebody_group.cpp)
target_include_directories(
	threebody_group PRIVATE ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(
	threebody_group PRIVATE event ROOT::Graf
)

add_executable(single_ppac_track_parameter single_ppac_track_parameter.cpp)
target_link_libraries(
	single_ppac_track_parameter PRIVATE Ceres::ceres event ROOT::Graf
)

# ppac_group
add_executable(ppac_group ppac_group.cpp)
target_link_libraries(ppac_group PRIVATE event ROOT::Core ROOT::RIO ROOT::Graf)