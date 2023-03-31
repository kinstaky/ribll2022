#include "include/detector/t0d2.h"

namespace ribll {

T0d2::T0d2(unsigned int run, const std::string &tag)
: Dssd(run, "t0d2", tag) {
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
