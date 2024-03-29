#ifndef __CSI_H__
#define __CSI_H__

#include <string>

#include "include/detector/detector.h"
#include "include/event/csi_event.h"

namespace ribll {

class GroupCsi : public Detector {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] name detector name
	/// @param[in] tag trigger tag
	/// @param[in] size number of csi in group
	///
	GroupCsi(
		unsigned int run,
		const std::string &name,
		const std::string &tag,
		unsigned int size
	);


	/// @brief default destructor
	///
	virtual ~GroupCsi() = default;

protected:
	unsigned int csi_size_;
};

class CircularCsi : public GroupCsi {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] name detector name
	/// @param[in] tag trigger tag
	///
	CircularCsi(
		unsigned int run,
		const std::string &name,
		const std::string &tag
	);


	/// @brief default destructor
	///
	virtual ~CircularCsi() = default;


	/// @brief extract trigger with detector events
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int MatchTrigger(
		double window_left,
		double window_right
	) override;
};


class SquareCsi : public GroupCsi {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] name detector name
	/// @param[in] tag trigger tag
	///
	SquareCsi(
		unsigned int run,
		const std::string &name,
		const std::string &tag
	);


	/// @brief default destructor
	///
	virtual ~SquareCsi() = default;


	/// @brief extract trigger with detector events
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int MatchTrigger(
		double window_left,
		double window_right
	) override;
};

}		// namespace ribll


#endif 		// __CSI_H__