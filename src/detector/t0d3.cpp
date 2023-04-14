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


int T0d3::NormalizeSides(TChain *chain, bool iteration) {
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
	const DssdFundamentalEvent&
) const {
	return true;
}

// bool T0D3::NormalizeFrontEnergyCheck
// (
// 	const CorrelatedEvent &correlation,
// 	bool iteration
// )
// {
// 	if
// 	(
// 		correlation.front_energy[0] < 2000
// 		|| correlation.back_energy[0] < 2000
// 	)
// 	{
// 		return false;
// 	}

// 	if (!iteration)
// 	{
// 	}
// 	else
// 	{
// 		if
// 		(
// 			abs
// 			(
// 				NormalizeEnergy
// 				(
// 					0,
// 					correlation.front_strip[0],
// 					correlation.front_energy[0]
// 				)
// 				- correlation.back_energy[0]
// 			)
// 			> 1000
// 		)
// 		{
// 			return false;
// 		}
// 	}
// 	return true;
// }


// bool T0D3::NormalizeBackEnergyCheck
// (
// 	const CorrelatedEvent &correlation,
// 	bool iteration
// )
// {
// 	if
// 	(
// 		correlation.front_energy[0] < 2000
// 		|| correlation.back_energy[0] < 2000
// 	)
// 	{
// 		return false;
// 	}

// 	if (!iteration)
// 	{
// 	}
// 	else
// 	{
// 		if
// 		(
// 			abs(
// 				NormalizeEnergy(
// 					0,
// 					correlation.front_strip[0],
// 					correlation.front_energy[0]
// 				)
// 				- NormalizeEnergy(
// 					1,
// 					correlation.back_strip[0],
// 					correlation.back_energy[0]
// 				)
// 			)
// 			> 1000
// 		)
// 		{
// 			return false;
// 		}
// 	}
// 	return true;
// }
}		// namespace ribll
