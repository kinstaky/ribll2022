#include "include/calculator/csi_energy_calculator.h"

#include <TMath.h>

namespace ribll {

namespace elc {

CsiEnergyCalculator::CsiEnergyCalculator(const std::string &projectile)
: projectile_(projectile)
, si_calculator_(projectile, "Si")
, al_calculator_(projectile, "Al")
, mylar_calculator_(projectile, "Mylar") {

}


double CsiEnergyCalculator::Energy(
	double theta,
	double si_energy,
	double si_thick,
	double dead_si_thick,
	double dead_al_thick,
	double mylar_thick
) const {
	// epsilon, stop if reach this standard
	double eps = 1e-4;
	// sensitive thickness in Si
	double si_sensitive_thick = si_thick - dead_si_thick * 2;
	// cos(theta)
	double cos_theta = TMath::Cos(theta);

	// lower bound of the total energy
	double lower_bound = si_calculator_.Energy(si_thick / cos_theta);
	// something went wrong if the Si energy is larger than the lower bound
	if (si_energy > lower_bound) return -1e5;

	// current distance of upper bound and lower bound
	double current_step = 1.0;
	// current hypothetical total energy, as the uppper bound
	double upper_bound;
	// search for lower and upper bound of total energy
	while (true) {
		upper_bound = lower_bound + current_step;
		// limited to the max energy
		if (upper_bound > si_calculator_.MaxEnergy()) {
			upper_bound = si_calculator_.MaxEnergy();
		}
		// incident range in Si under current total energy
		double si_range = si_calculator_.Range(upper_bound);
		// residual range after the sensitive layer
		si_range -= si_sensitive_thick / cos_theta;
		// residual energy after sensitive layer
		double residual_energy = si_calculator_.Energy(si_range);
		// loss energy under current total energy
		double loss_energy = upper_bound - residual_energy;
		if (loss_energy > si_energy) {
			// current total energy is too small
			// increase the lower bound and uppper bound
			lower_bound = upper_bound;
			// reach the maximum energy, the calculator can't handle
			if (lower_bound == si_calculator_.MaxEnergy()) return -2e5;
			current_step *= 2.0;
		} else {
			// get the lower and upper bound
			break;
		}
	}

	// max iteration to avoid infinite loop
	const int max_iteration = 1000;
	int current_iteration = 0;

	// current energy loss in sensitive layer
	double current_si_energy = 1000.0;
	// binary search for the total energy between the lower and upper bound
	while (fabs(current_si_energy - si_energy) > eps) {
		double current_total_energy = (upper_bound + lower_bound) / 2.0;
		// incident range in Si under current total energy
		double si_range = si_calculator_.Range(current_total_energy);
		// residual range after the sensitive layer
		si_range -= si_sensitive_thick / cos_theta;
		// residual energy after sensitive layer
		double residual_energy = si_calculator_.Energy(si_range);
		// loss energy undern current total energy
		current_si_energy = current_total_energy - residual_energy;
		// change the upper or lower bound
		if (current_si_energy > si_energy) {
			// current Si energy loss too large, increase lower bound
			lower_bound = current_total_energy;
		} else {
			// current Si energy loss too large, decrease upper bound
			upper_bound = current_total_energy;
		}
		// avoid infinite loop
		++current_iteration;
		if (current_iteration >= max_iteration) return -3e5;
	}

	// get the total energy
	double energy = (upper_bound + lower_bound) / 2.0;
	// incident range in Si under total energy
	double range = si_calculator_.Range(energy);
	// residual range after the sensitive layer and dead Si layer
	range -= (si_sensitive_thick + dead_si_thick) / cos_theta;
	// residual energy after sensitive layer and dead Si layer
	energy = si_calculator_.Energy(range);

	// incident range in Al material
	range = al_calculator_.Range(energy);
	// residual range after the dead Al layer
	range -= dead_al_thick / cos_theta;
	// residual energy after the dead Al layer
	energy = al_calculator_.Energy(range);

	// incident range in Mylar material
	range = mylar_calculator_.Range(energy);
	// residual range after the Mylar layer
	range -= mylar_thick / cos_theta;
	// residual energy after the Mylar layer
	energy = mylar_calculator_.Energy(range);

	return energy;
}

}	// namespace elc

}	// namespace ribll