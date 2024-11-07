# spectrum
add_executable(spectrum spectrum.cpp)
target_link_libraries(spectrum PRIVATE event ROOT::Hist d2_energy_calculator)

# spectrum V2
add_executable(spectrum_v2 spectrum_v2.cpp)
target_link_libraries(spectrum_v2 PRIVATE event ROOT::Hist d2_energy_calculator target_energy_calculator)

add_executable(spectrum_v3 spectrum_v3.cpp)
target_link_libraries(
	spectrum_v3 PRIVATE
	channel_v2_event target_energy_calculator csi_energy_calculator
	ROOT::Hist ROOT::RIO ROOT::Graf
)

# fit spectrum
add_executable(fit_spectrum fit_spectrum.cpp)
target_link_libraries(
	fit_spectrum
	PRIVATE
	ROOT::Hist
	event d2_energy_calculator target_energy_calculator
)

# roofit spectrum
add_executable(roofit_spectrum roofit_spectrum.cpp)
target_link_libraries(
	roofit_spectrum
	PRIVATE
	ROOT::Hist ROOT::RooFit ROOT::Gpad
	event d2_energy_calculator target_energy_calculator
	asymmetric_voigt background_poly
)

# spectrum 3a
add_executable(spectrum_3a spectrum_3a.cpp)
target_link_libraries(spectrum_3a PRIVATE event ROOT::Hist)

# spectrum 2a
add_executable(spectrum_2a spectrum_2a.cpp)
target_link_libraries(spectrum_2a PRIVATE event ROOT::Hist)

# spectrum 9Be
add_executable(spectrum_be9 spectrum_be9.cpp)
target_link_libraries(spectrum_be9 PRIVATE event ROOT::Hist)

# TAF edge spectrum
add_executable(taf_edge_spectrum taf_edge_spectrum.cpp)
target_link_libraries(taf_edge_spectrum PRIVATE event csi_energy_calculator ROOT::RIO)

# image recoil Q spectrum
add_executable(image_recoil_q image_recoil_q.cpp)
target_link_libraries(image_recoil_q PRIVATE channel_v2_event csi_energy_calculator)