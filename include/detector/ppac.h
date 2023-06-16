#ifndef __PPAC_H__
#define __PPAC_H__

#include <map>
#include <vector>

#include "include/detector/detector.h"
#include "include/event/ppac_event.h"


namespace ribll {

class Ppac : public Detector {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] name ppac name, xppac or vppac
	/// @param[in] tag trigger tag
	///
	Ppac(
		unsigned int run,
		const std::string &name,
		const std::string &tag
	);


	/// @brief match xia main trigger and build events
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int MatchTrigger(
		double window_left,
		double window_right
	) override;


	/// @brief check the sum of time and get the time difference
	/// @param double ignore
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Merge(double) override;


	/// @brief track PPAC events and calculate the reaction point
	/// @returns 0 if success, -1 otherwise
	virtual int Track();

private:

	int GetSumRange(PpacFundamentalEvent &fundamental, TTree *ipt, double *range);
};

}	// namespace ribll

#endif	// __PPAC_H__