#include "include/detector/t0d1.h"

namespace ribll {

// center of t0d1, in mm
const ROOT::Math::XYZVector t0d1_center{0.0, 0.0, 100.0};
// x range of t0d1, in mm
const std::pair<double, double> t0d1_x_range{-32.0, 32.0};
// y range of t0d1, in mm
const std::pair<double, double> t0d1_y_range{-32.0, 32.0};

T0d1::T0d1(unsigned int run, const std::string &tag)
: Dssd(run, "t0d1", tag) {

	center_ = t0d1_center;
	x_range_ = t0d1_x_range;
	y_range_ = t0d1_y_range;
}

//-----------------------------------------------------------------------------
//								geometry
//-----------------------------------------------------------------------------

ROOT::Math::XYZVector T0d1::CalculatePosition(double fs, double bs) const {
	double x = (x_range_.second - x_range_.first) / BackStrip();
	x = x * (bs + 0.5) + x_range_.first;
	double y = (y_range_.second - y_range_.first) / FrontStrip();
	y = y * (fs + 0.5) + y_range_.first;
	ROOT::Math::XYZVector result(x, y, 0.0);
	result += center_;
	return result;
}


//-----------------------------------------------------------------------------
//								normalize
//-----------------------------------------------------------------------------

int T0d1::NormalizeSides(TChain *chain, int iteration) {
	if (SideNormalize(chain, 0, 26, iteration)) {
		std::cerr << "Error: Normalize first side failed.\n";
		return -1;
	}
	if (SideNormalize(chain, 1, 35, iteration)) {
		std::cerr << "Error: Normalize second side failed.\n";
		return -1;
	}
	return 0;
}


bool T0d1::NormEnergyCheck(
	size_t,
	const DssdNormalizeEvent&
) const {
	return true;
}

}		// namespace ribll

