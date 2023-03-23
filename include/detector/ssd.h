#ifndef __SSD_H__
#define __SSD_H__

#include "include/detector/detector.h"

namespace ribll {

class Ssd : public Detector {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] name detector name
	///
	Ssd(unsigned int run, const std::string &name);


	/// @brief defautl destructor
	///
	virtual ~Ssd() = default;


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

}		// namespace ribll

#endif 		// __SSD_H__