#ifndef __DELTA_ENERGY_CALCULATOR_H__
#define __DELTA_ENERGY_CALCULATOR_H__

#include <string>
#include <vector>

#include <TFile.h>
#include <TSpline.h>

namespace ribll{

namespace elc {

class DeltaEnergyCalculator {
public:


	/// @brief constructor
	/// @param[in] telescope telescope name, e.g. t0, taf0
	/// @param[in] projectile projectile name, e.g. 1H, 4He
	///
	DeltaEnergyCalculator(
		const std::string &telescope,
		const std::string &projectile
	);


	/// @brief destructor
	///
	~DeltaEnergyCalculator();


	/// @brief initialize the delta energy calculator
	/// @param[in] telescope telescope name, e.g. t0
	/// @param[in] thickness thickness of Si in telescope, in um
	/// @param[in] projectiles name of projectiles, e.g. 1H, 4He 
	/// @returns 0 if success, -1 otherwise
	///
	static int Initialize(
		const std::string &telescope,
		const std::vector<double> &thickness,
		const std::vector<std::string> &projectiles
	);


	/// @brief calculate the energy loss in Si
	/// @param[in] layer layer to calculate, e.g. 0 for d1d2, 1 for d2d3...
	/// @param[in] delta_energy delta energy in the first layer in MeV
	///		(the first layer of d2d3 is d2)
	/// @returns energy loss in the second layer in MeV
	/// 	(the second layer of d2d3 is d3)
	///
	inline double Energy(unsigned short layer, double delta_energy) const {
		return e_de_funcs_[layer]->Eval(delta_energy);
	}


	/// @brief calculate the energy loss in Si
	/// @param[in] layer layer to calculate, e.g. 0 for d1d2, 1 for d2d3...
	/// @param[in] energy energy in the second layer in MeV
	///		(the second layer of d2d3 is d3)
	/// @returns energy in the first layer in MeV
	///		(the first layer of d2d3 is d2)
	///
	inline double DeltaEnergy(unsigned short layer, double energy) const {
		return de_e_funcs_[layer]->Eval(energy);
	}

private:
	// telescope name
	std::string telescope_;
	// projectile name
	std::string projectile_;
	// input file storing functions
	TFile *input_file_;
	// dE-E functions, evaluate dE from E
	std::vector<TSpline3*> de_e_funcs_;
	// E-dE functions, evaluate E from dE
	std::vector<TSpline3*> e_de_funcs_;
};

}	// namespace elc

}	// namespace ribll

#endif // __DELTA_ENERGY_CALCULATOR_H__