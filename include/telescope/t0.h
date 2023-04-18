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
	/// @param[in] angle_tolerance angle tolerance in tracking
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Track(double angle_tolerance) override;


	/// @brief identify particle in telescope
	/// @returns 0 if success, -1 otherwise
	///
	virtual int ParticleIdentify() override;


	/// @brief calibrate this telescope
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Calibrate() override;


	/// @brief get layers of Si detector
	/// @returns layers of Si detector
	///
	virtual inline size_t Layers() const {
		return 6;
	}


	/// @brief rebuild particle from layers of detectors
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Particle() override;

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
		const elc::CsiEnergyCalculator &csi_calculator
	) const;


	/// @brief calculate calibrated energy
	/// @param[in] layer 0-d1, 1-d2, 2-d3, 3-s1, 4-s2, 5-s3
	/// @param[in] energy energy in channel
	/// @returns calibrated energy in MeV
	///
	inline double CaliEnergy(size_t layer, double energy) const {
		return cali_params_[layer*2] + cali_params_[layer*2+1] * energy;
	}

};

};

#endif 	// __T0_H__