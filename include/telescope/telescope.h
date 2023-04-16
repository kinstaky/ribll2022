#ifndef __TELECOPE_H__
#define __TELECOPE_H__

#include <iostream>
#include <string>
#include <memory>

#include <TCutG.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "include/defs.h"
#include "include/energy_loss.h"

namespace ribll {

struct ParticleCut {
	unsigned short charge;
	unsigned short mass;
	std::unique_ptr<TCutG> cut;
};


class Telescope {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] name telescope name
	/// @param[in] tag trigger tag
	///
	Telescope(
		unsigned int run,
		const std::string &name,
		const std::string &tag
	);


	/// @brief default destructor
	///
	~Telescope() = default;


	/// @brief track particle in telescope
	/// @param[in] angle_tolerance angle tolerance in tracking
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Track(double angle_tolerance);


	/// @brief identify particle in telescope
	/// @returns 0 if success, -1 otherwise
	///
	virtual int ParticleIdentify();


	/// @brief calibrate this telescope
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Calibrate();


	/// @brief calibrate with alpha source
	/// @returns 0 if success, -1 otherwise
	///
	virtual int AlphaCalibrate();


	/// @brief calibrate csi with pid
	/// @returns 0 if success, -1 otherwise
	///
	virtual int CsiCalibrate();

protected:


	/// @brief calculate the energy in CsI based on the energy in last Si
	/// @param[in] projectile projectile particle, e.g. 1H, 4He
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
	double CalculateCsiEnergy(
		const std::string &projectile,
		double theta,
		double si_energy,
		double si_thick,
		double dead_si_thick = 0.5,
		double dead_al_thick = 0.3,
		double mylar_thick = 2.0
	);


	/// @brief read cut from file
	/// @param[in] prefix prefix of the cut
	/// @param[in] particle particle type
	/// @returns pointer to the cut if success, nullptr otherwise
	///
	std::unique_ptr<TCutG> ReadCut(
		const char *prefix,
		const char *particle
	) const;


	unsigned int run_;
	std::string name_;
	std::string tag_;
};

}		// namespace ribll

#endif		// __TELECOPE_H__