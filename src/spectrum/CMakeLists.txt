# asymmetric voigt
add_library(asymmetric_voigt asymmetric_voigt.cpp)
target_include_directories(
	asymmetric_voigt
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(asymmetric_voigt PUBLIC ROOT::MathCore)

# background polys
add_library(background_poly background_poly.cpp)
target_include_directories(
	background_poly
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(background_poly PUBLIC ROOT::MathCore)