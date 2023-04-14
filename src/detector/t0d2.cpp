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


int T0d2::NormalizeSides(TChain *chain, bool iteration) {
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
	const DssdFundamentalEvent&
) const {
	return true;
}

// bool T0D2::NormalizeFrontEnergyCheck(
// 	const CorrelatedEvent &correlation,
// 	bool iteration
// )
// {
// 	if
// 	(
// 		correlation.front_energy[0] < 5000
// 		|| correlation.back_energy[0] < 5000
// 	)
// 	{
// 		return false;
// 	}

// 	if (!iteration)
// 	{
// 		if
// 		(
// 			(
// 				correlation.front_energy[0] > 20000
// 				&& correlation.front_energy[0] < 34000
// 			)
// 			||
// 			(
// 				correlation.back_energy[0] > 22000
// 				&& correlation.back_energy[0] < 38000
// 			)
// 		)
// 		{
// 			return false;
// 		}
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


// bool T0D2::NormalizeBackEnergyCheck(
// 	const CorrelatedEvent &correlation,
// 	bool iteration
// )
// {
// 	if
// 	(
// 		correlation.front_energy[0] < 5000
// 		|| correlation.back_energy[0] < 5000
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

