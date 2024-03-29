#include "include/detector/csi.h"

namespace ribll {

GroupCsi::GroupCsi(
	unsigned int run,
	const std::string &name,
	const std::string &tag,
	unsigned int size
)
: Detector(run, name, tag)
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
	MatchTriggerStatistics &statistics
) {
	// initialize fundamental event
	fundamental_event.match = true;
	for (unsigned int i = 0; i < csi_size; ++i) {
		fundamental_event.time[i] = -1e5;
		fundamental_event.energy[i] = -1e5;
	}
	fundamental_event.cfd_flag = 0;

	// if any channel is valid
	bool valid = false;
	// used events
	unsigned int used_events = 0;
	// range of events time is equal to trigger time
	auto range = match_map.equal_range(trigger_time);
	for (auto iter = range.first; iter != range.second; ++iter) {
		if (fundamental_event.time[iter->second.index] > -9e4) {
			// Found more than one signal in one event in CsI(Tl).
			// So this is conflict event.
			fundamental_event.match = false;
			++statistics.conflict_events;
			break;
		} else {
			size_t index = iter->second.index;
			// Found the first signal in one CsI(Tl).
			fundamental_event.cfd_flag |=
				iter->second.cfd_flag ? (1 << index) : 0;
			fundamental_event.time[index] = iter->second.time - trigger_time;
			fundamental_event.energy[index] = iter->second.energy;
			fundamental_event.decode_entry[index] = iter->second.decode_entry;
			valid = true;
			++used_events;
		}
	}
	if (fundamental_event.match && valid) {
		++statistics.match_events;
		statistics.used_events += used_events;
	}
	if (!valid) fundamental_event.match = false;
	return;
}


CircularCsi::CircularCsi(
	unsigned int run,
	const std::string &name,
	const std::string &tag
)
: GroupCsi(run, name, tag, 12) {
}


int CircularCsi::MatchTrigger(
	double window_left,
	double window_right
) {
	return Detector::MatchTrigger<CsiMapEvent, CircularCsiFundamentalEvent>(
		window_left,
		window_right,
		FillEvent<CircularCsiFundamentalEvent, 12>
	);
}


SquareCsi::SquareCsi(
	unsigned int run,
	const std::string &name,
	const std::string &tag
)
: GroupCsi(run, name, tag, 4) {
}


int SquareCsi::MatchTrigger(
	double window_left,
	double window_right
) {
	return Detector::MatchTrigger<CsiMapEvent, SquareCsiFundamentalEvent>(
		window_left,
		window_right,
		FillEvent<SquareCsiFundamentalEvent, 4>
	);
}

}