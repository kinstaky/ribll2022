set(glog_DIR /data/d1/pwl/glog/lib64/cmake/glog/)
find_package(Ceres REQUIRED PATHS /data/d1/pwl/ceres-solver/)

# include directory
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/include")

add_executable(threebody_opt threebody_optimize.cpp ${PROJECT_INCLUDE_DIR}/optimize_utilities.h)
target_include_directories(
	threebody_opt PRIVATE ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(
	threebody_opt PRIVATE
	event delta_energy_calculator Ceres::ceres
)