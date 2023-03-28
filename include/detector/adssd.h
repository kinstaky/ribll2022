#ifndef __ADSSD_H__
#define __ADSSD_H__

#include <string>

#include "include/detector/dssd.h"
#include "include/event/dssd_event.h"

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


	//-------------------------------------------------------------------------
	//								gemometry
	//-------------------------------------------------------------------------

	/// @brief get front strip number
	/// @returns front strip number
	///
	virtual inline size_t FrontStrip() const {
		return 16;
	}


	/// @brief get back strip number
	/// @returns back strip number
	///
	virtual inline size_t BackStrip() const {
		return 8;
	}


	//-------------------------------------------------------------------------
	//							match trigger
	//-------------------------------------------------------------------------

	/// @brief match xia main trigger and build events
	/// @param[in] trigger_tag tag of trigger to chosse file
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int MatchTrigger(
		const std::string &trigger_tag,
		double window_left,
		double window_right
	) override;

};

}	// namespace ribll

#endif 		// __ADSSD_H__