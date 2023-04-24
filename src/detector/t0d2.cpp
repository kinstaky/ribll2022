#include "include/detector/t0d2.h"

namespace ribll {

// center of t0d2, in mm
const ROOT::Math::XYZVector t0d2_center{0.0, 0.0, 111.76};
// x range of t0d2, in mm
const std::pair<double, double> t0d2_x_range{-32.0, 32.0};
// y range of t0d2, in mm
const std::pair<double, double> t0d2_y_range{-32.0, 32.0};


T0d2::T0d2(unsigned int run, const std::string &tag)
: Dssd(run, "t0d2", tag) {

	center_ = t0d2_center;
	x_range_ = t0d2_x_range;
	y_range_ = t0d2_y_range;
}


int T0d2::NormalizeSides(TChain *chain, int iteration) {
	if (SideNormalize(chain, 0, 19, iteration)) {
		std::cerr << "Error: Normalize first side failed.\n";
		return -1;
	}
	if (SideNormalize(chain, 1, 14, iteration)) {
		std::cerr << "Error: Normalize second side failed.\n";
		return -1;
	}
	return 0;
}


bool T0d2::NormEnergyCheck(
	size_t,
	const DssdNormalizeEvent&
) const {
	return true;
}


}		// namespace ribll

