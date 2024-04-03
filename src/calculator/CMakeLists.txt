# include directory
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(CALCULATOR_INCLUDE_DIR ${PROJECT_INCLUDE_DIR}/calculator)

# range energy calculator library
add_library(
	range_energy_calculator
	range_energy_calculator.cpp
	${CALCULATOR_INCLUDE_DIR}/range_energy_calculator.h
)
target_include_directories(
	range_energy_calculator
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(
	range_energy_calculator
	PUBLIC ROOT::Core ROOT::RIO ROOT::Hist
)

# delta energy calculator library
add_library(
	delta_energy_calculator
	delta_energy_calculator.cpp
	${CALCULATOR_INCLUDE_DIR}/delta_energy_calculator.h
)
target_link_libraries(
	delta_energy_calculator
	PUBLIC range_energy_calculator
)

# D2 energy calculator library
add_library(
	d2_energy_calculator
	d2_energy_calculator.cpp
	${CALCULATOR_INCLUDE_DIR}/d2_energy_calculator.h
)
target_link_libraries(d2_energy_calculator PUBLIC range_energy_calculator)

# csi energy calculator library
add_library(
	csi_energy_calculator
	csi_energy_calculator.cpp
	${CALCULATOR_INCLUDE_DIR}/csi_energy_calculator.h
)
target_link_libraries(
	csi_energy_calculator
	PUBLIC range_energy_calculator
)

# target energy calculator library
add_library(
	target_energy_calculator
	target_energy_calculator.cpp
	${CALCULATOR_INCLUDE_DIR}/target_energy_calculator.h
)
target_link_libraries(
	target_energy_calculator
	PUBLIC range_energy_calculator
)