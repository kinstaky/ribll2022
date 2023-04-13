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

int T0d1::NormalizeSides(TChain *chain, bool iteration) {
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
	const DssdFundamentalEvent&
) const {
	return true;
}



// bool T0D1::NormalizeFrontEnergyCheck(
// 	const CorrelatedEvent &correlation,
// 	bool iteration
// )
// {
// 	if
// 	(
// 		correlation.front_energy[0] < 5000
// 		|| correlation.back_energy[0] < 5000
// 	)
// 		return false;


// 	if (!iteration)
// 	{
// 		if (correlation.front_strip[0] < 16)
// 		{
// 		}
// 		else if (correlation.front_strip[0] <= 32)
// 		{
// 			if
// 			(
// 				(
// 					correlation.front_energy[0] > 18000
// 					&& correlation.front_energy[0] < 30000
// 				)
// 				||
// 				(
// 					correlation.back_energy[0] > 30000
// 					&& correlation.back_energy[0] < 44000
// 				)
// 			)
// 			{
// 				return false;
// 			}

// 		}
// 		else if (correlation.front_strip[0] < 48)
// 		{
// 			if
// 			(
// 				(
// 					correlation.front_energy[0] > 23000
// 					&& correlation.front_energy[0] < 27000
// 				)
// 				||
// 				(
// 					correlation.back_energy[0] > 34000
// 					&& correlation.back_energy[0] < 41000
// 				)
// 			)
// 			{
// 				return false;
// 			}

// 			// if (
// 			// 	abs(correlation.front_energy[0] - norm_back_energy)
// 			// 	> 0.6 * norm_back_energy
// 			// ) return false;

// 		}
// 		else
// 		{

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
// 		) return false;

// 		if (correlation.front_strip[0] < 16)
// 		{


// 		}
// 		else if (correlation.front_strip[0] <= 32)
// 		{
// 			if
// 			(
// 				(
// 					correlation.front_energy[0] > 23000
// 					&& correlation.front_energy[0] < 26000
// 				)
// 				||
// 				(
// 					correlation.back_energy[0] > 34000
// 					&& correlation.back_energy[0] < 40000
// 				)
// 			)
// 			{
// 				return false;
// 			}


// 		}
// 		else if (correlation.front_strip[0] < 48)
// 		{
// 			if
// 			(
// 				(
// 					correlation.front_energy[0] > 17000
// 					&& correlation.front_energy[0] < 35000
// 				)
// 				||
// 				(
// 					correlation.back_energy[0] > 25000
// 					&& correlation.back_energy[0] < 49000
// 				)
// 			)
// 			{
// 				return false;
// 			}

// 		}
// 		else
// 		{
// 		}
// 	}

// 	return true;
// }



// bool T0D1::NormalizeBackEnergyCheck(const CorrelatedEvent &correlation, bool iteration) {
// 	if
// 	(
// 		correlation.front_energy[0] < 5000
// 		|| correlation.back_energy[0] < 5000
// 	)
// 		return false;

// 	if (!iteration)
// 	{
// 		// if (correlation.back_strip[0] < 16) {
// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 20000 && correlation.back_energy[0] < 22000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 30000 && correlation.back_energy[0] > 34000
// 		// 		)
// 		// 	)) return false;
// 		// } else if (correlation.back_strip[0] < 32) {
// 		// 	if (
// 		// 		abs(correlation.front_energy[0] - correlation.back_energy[0])
// 		// 		> 0.3 * correlation.back_energy[0]
// 		// 	) return false;

// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 20000 && correlation.back_energy[0] < 24000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 30000 && correlation.back_energy[0] > 38000
// 		// 		)
// 		// 	)) return false;
// 		// } else if (correlation.back_strip[0] < 48) {
// 		// 	if (
// 		// 		abs(correlation.front_energy[0] - correlation.back_energy[0])
// 		// 		> 0.5 * correlation.back_energy[0]
// 		// 	) return false;

// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 19000 && correlation.back_energy[0] < 28000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 29000 && correlation.back_energy[0] > 42000
// 		// 		)
// 		// 	)) return false;

// 		// } else {

// 		// }

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
// 			return false;


// 		// if (correlation.back_strip[0] < 16) {
// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 20000 && correlation.back_energy[0] < 22000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 30000 && correlation.back_energy[0] > 34000
// 		// 		)
// 		// 	)) return false;
// 		// } else if (correlation.back_strip[0] < 32) {
// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 20000 && correlation.back_energy[0] < 24000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 30000 && correlation.back_energy[0] > 38000
// 		// 		)
// 		// 	)) return false;
// 		// } else if (correlation.back_strip[0] < 48) {
// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 19000 && correlation.back_energy[0] < 28000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 29000 && correlation.back_energy[0] > 42000
// 		// 		)
// 		// 	)) return false;

// 		// } else {

// 		// }

// 	}

// 	return true;
// }



}		// namespace ribll

