#include "include/detector/tof.h"

#include <iostream>
#include <map>

#include <TH1F.h>

#include "include/defs.h"
#include "include/event/tof_event.h"


namespace ribll {


Tof::Tof(unsigned int run)
: Detector(run, "tof") {
}


/// @brief convert map event to fundamental event
/// @param[in] trigger_time trigger time to match
/// @param[in] match_map map_events order by trigger time
/// @param[out] fundamental_event converted fundamental event
/// @param[inout] statistics information about statistics
///
void FillEvent(
	double trigger_time,
	const std::multimap<double, TofMapEvent> &match_map,
	TofFundamentalEvent &fundamental_event,
	Detector::MatchTriggerStatistics &statistics
) {
	fundamental_event.time[0] = -1.0;
	fundamental_event.time[1] = -1.0;

	// check match events number
	size_t match_count = match_map.count(trigger_time);
	if (match_count == 1 || match_count == 2) {
		auto range = match_map.equal_range(trigger_time);
		for (auto iter = range.first; iter != range.second; ++iter) {
			fundamental_event.time[iter->second.index] = iter->second.time;
		}
		++statistics.match_events;
		statistics.used_events += match_count;
	} else if (match_count > 2) {
		++statistics.oversize_events;
	}
}


int Tof::MatchTrigger(double window_left, double window_right) {
	return Detector::MatchTrigger<
		TofMapEvent,
		TofFundamentalEvent,
		MatchTriggerStatistics
	>(
		window_left,
		window_right,
		FillEvent
	);
}

}		// namespace ribll