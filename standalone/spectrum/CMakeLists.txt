# spectrum
add_executable(spectrum spectrum.cpp)
target_link_libraries(spectrum PRIVATE event ROOT::Hist d2_energy_calculator)

# spectrum
add_executable(spectrum_v2 spectrum_v2.cpp)
target_link_libraries(spectrum_v2 PRIVATE event ROOT::Hist d2_energy_calculator)


# spectrum 3a
add_executable(spectrum_3a spectrum_3a.cpp)
target_link_libraries(spectrum_3a PRIVATE event ROOT::Hist)

# spectrum 2a
add_executable(spectrum_2a spectrum_2a.cpp)
target_link_libraries(spectrum_2a PRIVATE event ROOT::Hist)

# spectrum 9Be
add_executable(spectrum_be9 spectrum_be9.cpp)
target_link_libraries(spectrum_be9 PRIVATE event ROOT::Hist)