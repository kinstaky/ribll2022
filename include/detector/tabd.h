#ifndef __TAB_H__
#define __TAB_H__

#include "include/detector/adssd.h"

namespace ribll {

class Tabd : public Adssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] index index of taf, 0 to 5
	/// @param[in] tag trigger tag
	///
	Tabd(unsigned int run, unsigned int index, const std::string &tag);


	/// @brief default destructor
	///
	virtual ~Tabd() = default;


	//-------------------------------------------------------------------------
	//							match trigger
	//-------------------------------------------------------------------------

	/// @brief match xia main trigger and build events
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int MatchTrigger(
		double window_left,
		double window_right
	) override;
};

}

#endif		// __TAB_H__