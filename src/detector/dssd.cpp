#include "include/detector/dssd.h"

#include "include/event/dssd_event.h"


namespace ribll {

//-----------------------------------------------------------------------------
//						Dssd::MatchTriggerStatistics
//-----------------------------------------------------------------------------

Dssd::MatchTriggerStatistics::MatchTriggerStatistics(long long total)
: total_events(total)
, match_events(0)
, oversize_events(0) {
}


std::ostream& operator<<(
	std::ostream &os,
	const Dssd::MatchTriggerStatistics &st
) {
	return os << "events match rate "
		<< st.match_events << " / " << st.total_events << "  "
		<< double(st.match_events) / double(st.total_events) << "\n"
		<< st.oversize_events << " / " << st.total_events << "  "
		<< double(st.oversize_events) / double(st.total_events);
}

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
/// @param[inout] statistics information about statistics
///
void FillEvent(
	double trigger_time,
	const std::multimap<double, DssdMapEvent> &match_map,
	DssdFundamentalEvent &fundamental_event,
	Dssd::MatchTriggerStatistics &statistics
) {
	// initialize output fundamental event
	fundamental_event.front_hit = 0;
	fundamental_event.back_hit = 0;

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
		fundamental_event.front_hit = 0;
		fundamental_event.back_hit = 0;
		++statistics.oversize_events;
	} else {
		fundamental_event.front_hit = front_events.size();
		fundamental_event.back_hit = back_events.size();
	}
	for (unsigned short i = 0; i < fundamental_event.front_hit; ++i) {
		fundamental_event.front_strip[i] = front_events[i].strip;
		fundamental_event.front_time[i] = front_events[i].time;
		fundamental_event.front_energy[i] = front_events[i].energy;
	}
	for (unsigned short i = 0; i < fundamental_event.back_hit; ++i) {
		fundamental_event.back_strip[i] = back_events[i].strip;
		fundamental_event.back_time[i] = back_events[i].time;
		fundamental_event.back_energy[i] = back_events[i].energy;
	}

	if (
		fundamental_event.front_hit > 0
		&& fundamental_event.back_hit > 0
	) {
		++statistics.match_events;
	}

}

int Dssd::MatchTrigger(double window_left, double window_right) {
	return Detector::MatchTrigger<
		DssdMapEvent,
		DssdFundamentalEvent,
		MatchTriggerStatistics
	>(
		window_left,
		window_right,
		FillEvent
	);
}

//-----------------------------------------------------------------------------
//								T0d1
//-----------------------------------------------------------------------------

T0d1::T0d1(unsigned int run)
: Dssd(run, "t0d1") {
}



//-----------------------------------------------------------------------------
//								T0d2
//-----------------------------------------------------------------------------

T0d2::T0d2(unsigned int run)
: Dssd(run, "t0d2") {
}


//-----------------------------------------------------------------------------
//								T0d3
//-----------------------------------------------------------------------------

T0d3::T0d3(unsigned int run)
: Dssd(run, "t0d3") {
}

}