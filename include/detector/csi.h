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
	/// @param[in] size number of csi in group
	///
	GroupCsi(unsigned int run, const std::string &name, unsigned int size);


	/// @brief default destructor
	///
	virtual ~GroupCsi() = default;


	class MatchTriggerStatistics {
	public:
		long long total_events;
		long long match_events;
		long long oversize_events;

		/// @brief constructor
		/// @param[in] total total events
		///
		MatchTriggerStatistics(long long total);


		/// @brief overloaded operator<< function, output statistics
		/// @param[in] os ostream
		/// @param[in] statistics output object
		/// @returns ostream
		///
		friend std::ostream& operator<<(
			std::ostream &os,
			const MatchTriggerStatistics &statisics
		);
	};

protected:
	unsigned int csi_size_;
};

class CircularCsi : public GroupCsi {
public:
	
	/// @brief constructor
	/// @param[in] run run number 
	/// @param[in] name detector name
	///
	CircularCsi(unsigned int run, const std::string &name);


	/// @brief default destructor
	///
	virtual ~CircularCsi() = default;


	/// @brief match xia main trigger and build events
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int MatchTrigger(double window_left, double window_right) override;
};


class SquareCsi : public GroupCsi {
public:

	/// @brief constructor
	/// @param[in] run run number 
	/// @param[in] name detector name
	///
	SquareCsi(unsigned int run, const std::string &name);


	/// @brief default destructor
	///
	virtual ~SquareCsi() = default;


	/// @brief match xia main trigger and build events
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int MatchTrigger(double window_left, double window_right) override;

};

}		// namespace ribll


#endif 		// __CSI_H__