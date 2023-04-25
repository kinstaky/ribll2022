#ifndef __TAF_H__
#define __TAF_H__

#include "include/detector/adssd.h"

namespace ribll {

class Tafd : public Adssd {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] index index of taf, 0 to 5
	/// @param[in] tag trigger tag
	///
	Tafd(unsigned int run, unsigned int index, const std::string &tag);


	/// @brief default destructor
	///
	virtual ~Tafd() = default;


	/// @brief match xia main trigger and build events
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int MatchTrigger(
		double window_left,
		double window_right
	) override;


	/// @brief extract trigger with detector events
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int ExtractTrigger(
		double window_left,
		double window_right
	) override;


	/// @brief calibrate with alpha source
	/// @returns 0 if success, -1 otherwise
	///
	virtual int Calibrate() override;

private:
	// detector index, 0 to 5
	int index_;
};

}		// namespace ribll

#endif		// __TAF_H__