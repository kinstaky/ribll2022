#ifndef __TARGET_ENERGY_CALCULATOR_H__
#define __TARGET_ENERGY_CALCULATOR_H__

#include <string>

#include "include/calculator/range_energy_calculator.h"

namespace ribll {

namespace elc {


class TargetEnergyCalculator {
public:

	/// @brief constructor
	/// @param[in] projectile incident particle, e.g. 1H, 4He
	/// @param[in] material target material, e.g. C, CD2
	/// @param[in] density density of target, in mg/cm^2
	///
	TargetEnergyCalculator(
		const std::string &projectile,
		const std::string &material,
		double density
	);


	/// @brief default destructor
	virtual ~TargetEnergyCalculator() = default;

	double Energy(double depth, double energy);

private:
	RangeEnergyCalculator calculator_;
	double density_;
};

}	// elc

}	// ribll

#endif // __TARGET_ENERGY_CALCULATOR_H__