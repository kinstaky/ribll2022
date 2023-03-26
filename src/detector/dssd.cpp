#include "include/detector/dssd.h"

#include "include/event/dssd_event.h"


namespace ribll {

//-----------------------------------------------------------------------------
//								Dssd
//-----------------------------------------------------------------------------
Dssd::Dssd(unsigned int run, const std::string &name)
: Detector(run, name) {
}


/// @brief convert map event to fundamental event
/// @param[in] trigger_time trigger time to match
/// @param[in] match_map map_events order by trigger time
/// @param[out] fundamental_event converted fundamental event
/// @param[out] statistics information about statistics
/// @returns 0 if filling event is valid, -1 otherwise
///
int FillEvent(
	double trigger_time,
	const std::multimap<double, DssdMapEvent> &match_map,
	DssdFundamentalEvent &fundamental_event,
	MatchTriggerStatistics &statistics
) {
	// initialize output fundamental event
	fundamental_event.front_hit = 0;
	fundamental_event.back_hit = 0;
	fundamental_event.cfd_flag = 0;

	std::vector<DssdMapEvent> front_events;
	std::vector<DssdMapEvent> back_events;

	// check match events number
	auto range = match_map.equal_range(trigger_time);
	for (auto iter = range.first; iter != range.second; ++iter) {
		if (iter->second.side == 0) {
			front_events.push_back(iter->second);
		} else {
			back_events.push_back(iter->second);
		}
	}

	// sort events by strip
	std::sort(
		front_events.begin(),
		front_events.end(),
		[](const DssdMapEvent &x, const DssdMapEvent &y) {
			return x.strip < y.strip;
		}
	);
	std::sort(
		back_events.begin(),
		back_events.end(),
		[](const DssdMapEvent &x, const DssdMapEvent &y) {
			return x.strip < y.strip;
		}
	);

	// record events
	if (front_events.size() > 8 || back_events.size() > 8) {
		// Front or back events is more than 8, which is out of consideration.
		fundamental_event.front_hit = 0;
		fundamental_event.back_hit = 0;
		++statistics.oversize_events;
	} else {
		fundamental_event.front_hit = front_events.size();
		fundamental_event.back_hit = back_events.size();
	}

	// check multiple hit on the same strip
	bool strip_multiple_hit = false;
	for (size_t i = 1; i < front_events.size(); ++i) {
		if (front_events[i].strip == front_events[i-1].strip) {
			strip_multiple_hit = true;
			fundamental_event.front_hit = 0;
		}
	}
	for (size_t i = 1; i < back_events.size(); ++i) {
		if (back_events[i].strip == back_events[i-1].strip) {
			strip_multiple_hit = true;
			fundamental_event.back_hit = 0;
		}
	}
	statistics.conflict_events += strip_multiple_hit ? 1 : 0;


	for (unsigned short i = 0; i < fundamental_event.front_hit; ++i) {
		fundamental_event.front_strip[i] = front_events[i].strip;
		fundamental_event.front_time[i] = front_events[i].time - trigger_time;
		fundamental_event.front_energy[i] = front_events[i].energy;
		fundamental_event.cfd_flag |=
			front_events[i].cfd_flag ? (1 << i) : 0;
	}
	for (unsigned short i = 0; i < fundamental_event.back_hit; ++i) {
		fundamental_event.back_strip[i] = back_events[i].strip;
		fundamental_event.back_time[i] = back_events[i].time - trigger_time;
		fundamental_event.back_energy[i] = back_events[i].energy;
		fundamental_event.cfd_flag |=
			back_events[i].cfd_flag ? (1 << (i + 8)) : 0;
	}

	if (
		fundamental_event.front_hit > 0
		&& fundamental_event.back_hit > 0
	) {
		++statistics.match_events;
		statistics.used_events +=
			fundamental_event.front_hit + fundamental_event.back_hit;
		return 0;
	}

	return -1;
}


void FillEventInMatch(
	double trigger_time,
	const std::multimap<double, DssdMapEvent> &match_map,
	DssdFundamentalEvent &fundamental_event,
	MatchTriggerStatistics &statistics
) {
	FillEvent(trigger_time, match_map, fundamental_event, statistics);
}


int FillEventInExtract(
	double trigger_time,
	const std::multimap<double, DssdMapEvent> &match_map,
	DssdFundamentalEvent &fundamental_event,
	std::vector<MatchTriggerStatistics> &statistics
) {
	return
		FillEvent(trigger_time, match_map, fundamental_event, statistics[0]);
}


int Dssd::MatchTrigger(
	const std::string &trigger_tag,
	double window_left,
	double window_right
) {
	return Detector::MatchTrigger<DssdMapEvent, DssdFundamentalEvent>(
		trigger_tag,
		window_left,
		window_right,
		FillEventInMatch
	);
}


int Dssd::ExtractTrigger(
	const std::string &trigger_tag,
	double window_left,
	double window_right
) {
	return Detector::ExtractTrigger<DssdMapEvent, DssdFundamentalEvent>(
		trigger_tag,
		{name_},
		window_left,
		window_right,
		FillEventInExtract
	);
}

}