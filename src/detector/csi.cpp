#include "include/detector/csi.h"

namespace ribll {

GroupCsi::GroupCsi(
	unsigned int run,
	const std::string &name,
	unsigned int size
)
: Detector(run, name)
, csi_size_(size) {
}


/// @brief convert map event to fundamental event
/// @param[in] trigger_time trigger time to match
/// @param[in] match_map map_events order by trigger time
/// @param[out] fundamental_event converted fundamental event
/// @param[inout] statistics information about statistics
///
template<typename FundamentalEvent, unsigned int csi_size>
void FillEvent(
	double trigger_time,
	const std::multimap<double, CsiMapEvent> &match_map,
	FundamentalEvent &fundamental_event,
	Detector::MatchTriggerStatistics &statistics
) {
	// initialize fundamental event
	fundamental_event.match = true;
	for (unsigned int i = 0; i < csi_size; ++i) {
		fundamental_event.time[i] = -1.0;
		fundamental_event.energy[i] = -1.0;
	}
	// check match events number
	// if any channel is valid
	bool valid = false;
	auto range = match_map.equal_range(trigger_time);
	for (auto iter = range.first; iter != range.second; ++iter) {
		if (fundamental_event.time[iter->second.index] > 0.0) {
			fundamental_event.match = false;
			++statistics.oversize_events;
		} else {
			fundamental_event.time[iter->second.index] = iter->second.time;
			fundamental_event.energy[iter->second.index] = iter->second.energy;
			valid = true;
		}
	}
	if (fundamental_event.match && valid) {
		++statistics.match_events;
	}
	return;
}


CircularCsi::CircularCsi(unsigned int run, const std::string &name)
: GroupCsi(run, name, 12) {
}


int CircularCsi::MatchTrigger(double window_left, double window_right) {
	return Detector::MatchTrigger<
		CsiMapEvent,
		CircularCsiFundamentalEvent,
		MatchTriggerStatistics
	>(
		window_left,
		window_right,
		FillEvent<CircularCsiFundamentalEvent, 12>
	);
}


SquareCsi::SquareCsi(unsigned int run, const std::string &name)
: GroupCsi(run, name, 4) {
}


int SquareCsi::MatchTrigger(double window_left, double window_right) {
	return Detector::MatchTrigger<
		CsiMapEvent,
		SquareCsiFundamentalEvent,
		MatchTriggerStatistics
	>(
		window_left,
		window_right,
		FillEvent<SquareCsiFundamentalEvent, 4>
	);
}

}