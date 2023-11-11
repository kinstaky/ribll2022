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
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Track();


	/// @brief slice track
	/// @returns 0 if success, -1 otherwise
	///
	virtual int SliceTrack();


	/// @brief identify particle in telescope
	/// @returns 0 if success, -1 otherwise
	///
	virtual int ParticleIdentify();


	/// @brief calibrate this telescope
	/// @param[in] end_run end of run to chain, inclusive
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Calibrate(unsigned int end_run);


	/// @brief calibrate with alpha source
	/// @returns 0 if success, -1 otherwise
	///
	virtual int AlphaCalibrate();


	/// @brief calibrate csi with pid
	/// @param[in] end_run end of run to chain, inclusive
	/// @returns 0 if success, -1 otherwise
	///
	virtual int CsiCalibrate(unsigned int end_run);


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


	/// @brief generate calibrate result
	/// @returns 0 if success, -1 otherwise
	///
	virtual int CalibrateResult();

protected:

	/// @brief read cut from file
	/// @param[in] name cut name
	/// @returns pointer to the cut if success, nullptr otherwise
	///
	std::unique_ptr<TCutG> ReadCut(const char *name) const;

	// run number
	unsigned int run_;
	// telescope name, e.g. t0, taf0
	std::string name_;
	// trigger tag
	std::string tag_;
};


}		// namespace ribll

#endif		// __TELECOPE_H__