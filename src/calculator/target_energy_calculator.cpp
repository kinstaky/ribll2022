#include "include/calculator/target_energy_calculator.h"

namespace ribll {

namespace elc {

TargetEnergyCalculator::TargetEnergyCalculator(
	const std::string &projectile,
	const std::string &material,
	double density
)
: calculator_(projectile, material)
, density_(density)
{}

double TargetEnergyCalculator::Energy(double depth, double energy) {
	double range = calculator_.Range(energy);
	range -= depth * density_;
	return calculator_.Energy(range);
}

}	// elc

}	// ribll