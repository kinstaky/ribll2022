#include "include/detector/ssd.h"

#include "include/event/ssd_event.h"

namespace ribll {

Ssd::Ssd(unsigned int run, const std::string &name)
: Detector(run, name) {
}


/// @brief convert map event to fundamental event
/// @param[in] trigger_time trigger time to match
/// @param[in] match_map map_events order by trigger time
/// @param[out] fundamental_event converted fundamental event
/// @param[out] statistics information about statistics
/// @returns 0 if filling event is valid, -1 otherwise
///
void FillEvent(
	double trigger_time,
	const std::multimap<double, SsdEvent> &match_map,
	SsdEvent &fundamental_event,
	MatchTriggerStatistics &statistics
) {
	// initialize output fundamnetal event
	fundamental_event.cfd_flag = 0;
	fundamental_event.time = -1e5;
	fundamental_event.energy = 0;

	size_t match_count = match_map.count(trigger_time);
	if (match_count > 1) {
		++statistics.oversize_events;
	} else if (match_count == 1) {
		const auto &event = match_map.equal_range(trigger_time).first->second;
		fundamental_event.cfd_flag = event.cfd_flag;
		fundamental_event.time = event.time - trigger_time;
		fundamental_event.energy = event.energy;
		++statistics.match_events;
		++statistics.used_events;
	}
}


int Ssd::MatchTrigger(
	const std::string &trigger_tag,
	double window_left,
	double window_right
) {
	return Detector::MatchTrigger<SsdEvent, SsdEvent>(
		trigger_tag,
		window_left,
		window_right,
		FillEvent
	);
}


}		// namespace ribll