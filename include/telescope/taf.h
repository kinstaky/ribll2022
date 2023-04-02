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
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Track() override;


	/// @brief calibrate this telescope
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Calibrate() override;


	/// @brief calibrate with alpha source
	/// @returns 0 if success, -1 otherwise
	///
	virtual int AlphaCalibrate() override;

private:
	unsigned int index_;
};

}		// namespace ribll

#endif		// __TAF_H__