#ifndef __T0_H__
#define __T0_H__

#include "include/telescope/telescope.h"

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
private:

};

};

#endif 	// __T0_H__