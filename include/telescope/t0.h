#ifndef __T0_H__
#define __T0_H__

#include "include/telescope/telescope.h"
#include "include/event/t0_event.h"
#include "include/event/particle_type_event.h"
#include "include/event/particle_event.h"

namespace ribll {

class T0 : public Telescope {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] tag trigger tag
	///
	T0(unsigned int run, const std::string &tag);


	/// @brief track particle
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Track() override;


	/// @brief identify particle in telescope
	/// @returns 0 if success, -1 otherwise
	///
	virtual int ParticleIdentify() override;


	/// @brief calibrate this telescope
	/// @param[in] end_run end of run to chain, inclusive
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Calibrate(unsigned int end_run) override;


	/// @brief read calibrate parameters from file
	/// @returns 0 if success, -1 otherwise
	///
	int ReadCalibrateParameters(unsigned int run = 9999);


	/// @brief write calibrate parameters to file
	/// @returns 0 if success, -1 otherwise
	///
	int WriteCalibrateParameters(unsigned int run = 9999) const;


	/// @brief generate calibrate result
	/// @returns 0 if success, -1 otherwise
	///
	virtual void CalibrateResult(T0Event &t0_event);


	/// @brief show calibration result
	/// @returns 0 if success, -1 otherwise
	///
	virtual int ShowCalibration();


	/// @brief get layers of Si detector
	/// @returns layers of Si detector
	///
	virtual inline size_t Layers() const {
		return 6;
	}


	/// @brief rebuild particle from layers of detectors
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Rebuild() override;


	/// @brief merge and track at the same time
	/// @returns 0 if success, -1 otherwise
	///
	virtual int MergeAndTrack();

private:

	/// @brief calculate the total kinetic energy lost in telescope
	/// @param[in] t0 t0 telescope event
	/// @param[in] type particle type event
	/// @param[in] index index of the particle in events to calculate
	/// @param[in] csi_calculator CsI energy calculator
	/// @returns calculated total kinetic energy
	///
	double TotalEnergy(
		const T0Event &t0,
		const ParticleTypeEvent &type,
		const size_t index,
		const elc::CsiEnergyCalculator &csi_calculator,
		const elc::DeltaEnergyCalculator &delta_calculator
	) const;


	/// @brief calculate calibrated energy
	/// @param[in] layer 0-d1, 1-d2, 2-d3, 3-s1, 4-s2, 5-s3
	/// @param[in] energy energy in channel
	/// @returns calibrated energy in MeV
	///
	inline double CaliEnergy(size_t layer, double energy) const {
		return cali_params_[layer*2] + cali_params_[layer*2+1] * energy;
	}

	// calibrate parameters
	double cali_params_[12];
};

};

#endif 	// __T0_H__