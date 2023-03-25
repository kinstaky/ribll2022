#ifndef __TOF_H__
#define __TOF_H__

#include <vector>
#include <ostream>

#include <TFile.h>
#include <TTree.h>

#include "include/detector/detector.h"

namespace ribll {

/// class for ToF detectors T1 and T2.
class Tof : public Detector {
public:

	/// @brief constructor
	/// @param[in] run run number
	///
	Tof(unsigned int run);


	///	@brief default destructor
	///
	~Tof() = default;


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


	// /// @brief extract trigger with detector events
	// /// @param[in] trigger_tag extract from trigger with this tag
	// /// @param[in] window_left left edge of match window
	// /// @param[in] window_right right edge of match window
	// /// @returns 0 if success, -1 otherwise
	// ///
	// virtual int ExtractTrigger(
	// 	const std::string &trigger_tag,
	// 	double window_left,
	// 	double window_right
	// ) override;


	/// @brief identify beam type by ToF
	/// @returns 0 if success, -1 otherwise
	///
	virtual int BeamIdentify();
};


}	// namespace ribll

#endif		// __TOF_H__