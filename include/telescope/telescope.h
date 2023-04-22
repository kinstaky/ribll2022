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
#include "include/calculator/range_energy_calculator.h"
#include "include/calculator/delta_energy_calculator.h"
#include "include/calculator/csi_energy_calculator.h"
#include "include/event/particle_type_event.h"
#include "include/event/particle_event.h"

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
	/// @param[in] length number of runs to chain
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Calibrate(unsigned int length);


	/// @brief calibrate with alpha source
	/// @returns 0 if success, -1 otherwise
	///
	virtual int AlphaCalibrate();


	/// @brief calibrate csi with pid
	/// @param[in] length number of runs to chain
	/// @returns 0 if success, -1 otherwise
	///
	virtual int CsiCalibrate(unsigned int length);


	/// @brief get layers of Si detector
	/// @returns layers of Si detector
	///
	virtual inline size_t Layers() const {
		return 1;
	}


	/// @brief rebuild the particle from layers of detectors
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Rebuild();

protected:

	/// @brief read cut from file
	/// @param[in] prefix prefix of the cut
	/// @param[in] particle particle type
	/// @returns pointer to the cut if success, nullptr otherwise
	///
	std::unique_ptr<TCutG> ReadCut(
		const char *prefix,
		const char *particle
	) const;


	/// @brief read calibrate parameters from file
	/// @returns 0 if success, -1 otherwise
	///
	int ReadCalibrateParameters();


	/// @brief write calibrate parameters to file
	/// @returns 0 if success, -1 otherwise
	///
	int WriteCalibrateParameters() const;


	// run number
	unsigned int run_;
	// telescope name, e.g. t0, taf0
	std::string name_;
	// trigger tag
	std::string tag_;
	// calibrate parameters
	double cali_params_[12];
};


}		// namespace ribll

#endif		// __TELECOPE_H__