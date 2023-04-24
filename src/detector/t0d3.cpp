#include "include/detector/t0d3.h"

namespace ribll {

// center of t0d3, in mm
const ROOT::Math::XYZVector t0d3_center{0.0, 0.0, 123.52};
// x range of t0d3, in mm
const std::pair<double, double> t0d3_x_range{-32.0, 32.0};
// y range of t0d3, in mm
const std::pair<double, double> t0d3_y_range{-32.0, 32.0};


T0d3::T0d3(unsigned int run, const std::string &tag)
: Dssd(run, "t0d3", tag) {

	center_ = t0d3_center;
	x_range_ = t0d3_x_range;
	y_range_ = t0d3_y_range;
}


int T0d3::NormalizeSides(TChain *chain, int iteration) {
	if (SideNormalize(chain, 0, 20, iteration)) {
		std::cerr << "Error: Normalize first side failed.\n";
		return -1;
	}
	if (SideNormalize(chain, 1, 13, iteration)) {
		std::cerr << "Error: Normalize second side failed.\n";
		return -1;
	}
	return 0;
}


bool T0d3::NormEnergyCheck(
	size_t,
	const DssdNormalizeEvent&
) const {
	return true;
}

}		// namespace ribll
