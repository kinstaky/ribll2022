#ifndef __TELECOPE_H__
#define __TELECOPE_H__

#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <string>

#include "include/defs.h"

namespace ribll {

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

	/// @brief calibrate this telescope
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Calibrate();


	/// @brief calibrate with alpha source
	/// @returns 0 if success, -1 otherwise
	///
	virtual int AlphaCalibrate();

protected:
	unsigned int run_;
	std::string name_;
	std::string tag_;
};

}		// namespace ribll

#endif		// __TELECOPE_H__