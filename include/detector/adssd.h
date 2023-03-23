#ifndef __ADSSD_H__
#define __ADSSD_H__

#include <string>

#include "include/detector/dssd.h"

namespace ribll {

class Adssd : public Dssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] name detector name
	///
	Adssd(unsigned int run, const std::string &name);


	/// @brief default destructor
	///
	virtual ~Adssd() = default;
};


class Taf : public Adssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] index index of taf, 0 to 5
	///
	Taf(unsigned int run, unsigned int index);


	/// @brief default destructor
	///
	virtual ~Taf() = default;


	/// @brief extract trigger with detector events
	/// @param[in] trigger_tag extract from trigger with this tag
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int ExtractTrigger(
		const std::string &trigger_tag,
		double window_left,
		double window_right
	) override;
};


class Tab : public Adssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] index index of taf, 0 to 5
	///
	Tab(unsigned int run, unsigned int index);


	/// @brief default destructor
	///
	virtual ~Tab() = default;
};


}	// namespace ribll

#endif 		// __ADSSD_H__