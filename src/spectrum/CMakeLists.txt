# asymmetric voigt
add_library(asymmetric_voigt asymmetric_voigt.cpp)
target_include_directories(
	asymmetric_voigt
	PUBLIC ${PROJECT_SOURCE_DIR} ROOT_INCLUDE_DIRS
)
target_link_libraries(asymmetric_voigt PUBLIC ROOT::MathCore)