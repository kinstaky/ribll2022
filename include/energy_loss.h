#ifndef __ENERGY_LOSS_H__
#define __ENERGY_LOSS_H__

#include <string>

#include <TFile.h>
#include <TSpline.h>

namespace ribll {

// namespace energy loss
namespace el {

typedef std::pair<std::string, std::string> ProjectileMaterial;


/// @brief convert the projectile material function from LISE++ to TSpline3
/// @param[in] list list of projectile material combination
/// @returns 0 if success, -1 otherwise
///
int Initialize(const std::vector<ProjectileMaterial> &list);


class EnergyLossCalculator {
public:

	/// @brief constructor
	/// @param[in] projectile incident particle, e.g. 1H, 4He
	/// @param[in] material loss in material, e.g. Si, Al, Mylar
	///
	EnergyLossCalculator(
		const std::string &projectile,
		const std::string &material
	);


	/// @brief destructor
	///
	~EnergyLossCalculator();


	/// @brief get mass number of the projectile
	/// @returns mass number
	inline unsigned int Mass() const {
		return mass_;
	}


	/// @brief calculate the incident range (um) from energy (MeV)
	/// @param[in] energy incident energy in MeV
	/// @returns incident range in um
	///
	double Range(double energy) const;


	/// @brief calculate the incident energy (MeV) from range (um)
	/// @param[in] range incident range in um
	/// @returns incident energy in MeV
	///
	double Energy(double range) const;


	/// @brief get maximun energy in the calculator 
	/// @returns maximun energy in the calculator
	///
	inline double MaxEnergy() const {
		return mass_ * 500.0;
	}

private:
	// projectile particle
	std::string projectile_;
	// material
	std::string material_;
	// mass number
	unsigned int mass_;
	// input file storing functions
	TFile *input_file_;
	// energy-range function, evaluate energy from range
	TSpline3 *energy_range_func_;
	// range-energy function, evaluate range from energy
	TSpline3 *range_energy_func_;
};

}		// namespace el

}		// namespace ribll

#endif		// __ENERGY_LOSS_H__