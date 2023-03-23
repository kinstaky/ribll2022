#ifndef __VME_TRIGGER_H__
#define __VME_TRIGGER_H__

#include "include/detector/detector.h"

namespace ribll {

class VmeTrigger : public Detector {
public:

	/// @brief constructor
	/// @param[in] run run number
	///
	VmeTrigger(unsigned int run);


	/// @brief default destructor
	///
	virtual ~VmeTrigger() = default;


	/// @brief match xia main trigger and build events
	/// @param[in] trigger_tag tag of trigger to chosse file
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

#endif		// __VME_TRIGGER_H__