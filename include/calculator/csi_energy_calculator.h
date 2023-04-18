#ifndef __CSI_ENERGY_CALCULATOR_H__
#define __CSI_ENERGY_CALCULATOR_H__

#include "include/calculator/range_energy_calculator.h"

namespace ribll {

namespace elc {

class CsiEnergyCalculator {
public:
	/// @brief constructor
	/// @param[in] projectile incident particle, e.g. 1H, 4He
	///
	CsiEnergyCalculator(const std::string &projectile);


	/// @brief default destructor
	///
	~CsiEnergyCalculator() = default;


	/// @brief calculate the energy in CsI(Tl) based on the energy in last Si
	/// @param[in] theta theta of projectile particle
	/// @param[in] si_energy energy loss in last Si, in MeV
	/// @param[in] si_thick base thickness of last Si, in um
	/// @param[in] dead_si_thick dead Si thickness of last Si, in um
	/// @param[in] dead_al_thick dead Al thickness of last Si, in um
	/// @param[in] mylar_thick Mylar thickness of CsI, in um
	/// @returns energy should be lost in CsI if success,
	/// 	-1e5 for energy can't pierce the Si
	///		-2e5 for energy larger the maximum energy in calculator
	///		-3e5 for iteration over maximun iterations
	///
	double Energy(
		double theta,
		double si_energy,
		double si_thick,
		double dead_si_thick = 0.5,
		double dead_al_thick = 0.3,
		double mylar_thick = 2.0
	) const;

private:
	// projectile particle
	std::string projectile_;
	// range-energy calculator to calculate energy lost in Si
	RangeEnergyCalculator si_calculator_;
	// range-energy calculator to calculate energy lost in Al
	RangeEnergyCalculator al_calculator_;
	// range-energy calculator to calculate energy lost in Mylar
	RangeEnergyCalculator mylar_calculator_;
};

}	// namespace elc

}	// namespace ribll

#endif // __CSI_ENERGY_CALCULATOR_H__