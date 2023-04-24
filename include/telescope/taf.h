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
	/// @param[in] end_run end of run to chain, inclusive
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Calibrate(unsigned int end_run) override;


	/// @brief calibrate csi with pid
	/// @param[in] end_run end of run to chain, inclusive
	/// @returns 0 if success, -1 otherwise
	///
	virtual int CsiCalibrate(unsigned int end_run) override;


	/// @brief rebuild the particle from layers of detectors
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Rebuild() override;

protected:

	/// @brief read CsI calibrate parameters from file
	/// @returns 0 if success, -1 otherwise
	///
	int ReadCsiCalibrateParameters();


	/// @brief write CsI calibrate parameters to file
	/// @returns 0 if success, -1 otherwise
	///
	int WriteCsiCalibrateParameters();

private:
	unsigned int index_;
	double csi_calibrate_params_[4][4][3];
};

}		// namespace ribll

#endif		// __TAF_H__