#ifndef __D2_ENERGY_CALCULATOR_H__
#define __D2_ENERGY_CALCULATOR_H__

#include <TFile.h>
#include <TSpline.h>

#include "include/calculator/range_energy_calculator.h"

namespace ribll {

namespace elc {

class D2EnergyCalculator {
public:

	/// @brief constructor
	/// @param[in] projectile incident particle, e.g. 1H, 4He
	///
	D2EnergyCalculator(const std::string &projectile);


	/// @brief destructor
	///
	~D2EnergyCalculator();


	/// @brief initialize the delta energy calculator
	/// @param[in] thickness thickness of Si in telescope, in um
	/// @param[in] projectiles name of projectiles, e.g. 1H, 4He
	/// @returns 0 if success, -1 otherwise
	///
	static int Initialize(
		const std::vector<double> &thickness,
		const std::vector<std::string> &projectiles
	);


	/// @brief calculate the energy loss in Si
	/// @param[in] layer layer to calculate, e.g. 0 for d1d2, 1 for d2d3...
	/// @param[in] delta_energy delta energy in the first layer in MeV
	///		(the first layer of d2d3 is d2)
	/// @param[in] stop true if stop in the second layer
	/// @returns energy loss in the second layer in MeV
	/// 	(the second layer of d2d3 is d3)
	///
	double Energy(
		unsigned short layer,
		double delta_energy,
		double theta,
		bool stop = true
	) const;
	//  {
	// 	if (stop) return e_de_stop_funcs_[layer]->Eval(delta_energy);
	// 	return e_de_pass_funcs_[layer]->Eval(delta_energy);
	// }


	/// @brief calculate the energy loss in Si
	/// @param[in] layer layer to calculate, e.g. 0 for d1d2, 1 for d2d3...
	/// @param[in] energy energy in the first layer in MeV
	///		(the first layer of d2d3 is d2)
	/// @param[in] stop true if stop in the second layer
	/// @returns energy loss in the first layer in MeV
	/// 	(the first layer of d2d3 is d2)
	///
	double DeltaEnergy(
		unsigned short layer,
		double energy,
		double theta,
		bool stop = true
	) const;
	//  {
	// 	if (stop) return de_e_stop_funcs_[layer]->Eval(energy);
	// 	return de_e_pass_funcs_[layer]->Eval(energy);
	// }

private:
	std::vector<TSpline3*> de_e_stop_funcs_;
	std::vector<TSpline3*> de_e_pass_funcs_;
	std::vector<TSpline3*> e_de_stop_funcs_;
	std::vector<TSpline3*> e_de_pass_funcs_;
	TFile *input_file_;

	RangeEnergyCalculator si_calculator_;
};

}	// elc

}	// ribll

#endif	// __D2_ENERGY_CALCULATOR_H__