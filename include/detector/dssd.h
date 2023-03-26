#ifndef __DSSD_H__
#define __DSSD_H__

#include <string>

#include "include/detector/detector.h"

namespace ribll {

class Dssd : public Detector {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] name detector name
	///
	Dssd(unsigned int run, const std::string &name);


	/// @brief default destructor
	///
	virtual ~Dssd() = default;


	//-------------------------------------------------------------------------
	//								gemometry
	//-------------------------------------------------------------------------

	/// @brief get front strip number
	/// @returns front strip number
	///
	virtual inline size_t FrontStrip() const {
		return 32;
	}


	/// @brief get back strip number
	/// @returns back strip number
	///
	virtual inline size_t BackStrip() const {
		return 32;
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

}		// namespace ribll

#endif		// __DSSD_H__