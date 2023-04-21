#ifndef __TAF_H__
#define __TAF_H__

#include "include/telescope/telescope.h"

namespace ribll {

class Taf : public Telescope {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] index telescope index, 0~5
	/// @param[in] tag trigger tag
	///
	Taf(unsigned int run, unsigned int index, const std::string &tag);


	/// @brief track particle
	/// @param[in] angle_tolerance angle tolerance
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


	/// @brief calibrate csi with pid
	/// @returns 0 if success, -1 otherwise
	///
	virtual int CsiCalibrate() override;


	/// @brief rebuild the particle from layers of detectors
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Rebuild() override;

private:
	unsigned int index_;
};

}		// namespace ribll

#endif		// __TAF_H__