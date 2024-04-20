#ifndef __LOST_ENERGY_CALCULATOR_H__
#define __LOST_ENERGY_CALCULATOR_H__

#include <string>

#include <TFile.h>
#include <TSpline.h>

namespace ribll {

namespace elc {

class LostEnergyCalculator {
public:

	/// @brief constructor
	/// @param[in] projectile projectile name, e.g. 1H, 4He
	///
	LostEnergyCalculator(const std::string &projectile);

	
	/// @brief destructor
	///
	~LostEnergyCalculator();



	static int Initialize(const std::vector<std::string> &projectiles);


	/// @brief get energy lost (dE/dx) from energy (E)
	///
	inline double EnergyLost(double energy, double depth) {
		return dedx_from_e_func_->Eval(energy) * depth;
	}


private:
	// projectile name
	std::string projectile_;
	// input file
	TFile *input_file_;
	// energy lost from energy (dE/dx from E)
	TSpline3 *dedx_from_e_func_;
};

}	// elc

}	// ribll

#endif	// __LOST_ENERGY_CALCULATOR_H__