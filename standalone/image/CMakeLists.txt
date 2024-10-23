add_executable(taf_pid_image taf_pid_image.cpp)
target_link_libraries(taf_pid_image PRIVATE event ROOT::RIO ROOT::Tree ROOT::Graf ROOT::Gpad)