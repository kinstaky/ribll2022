#include "include/statistics/match_trigger_statistics.h"

namespace ribll {

MatchTriggerStatistics::MatchTriggerStatistics(
	unsigned int run,
	const std::string &detector,
	const std::string &tag,
	const std::string &extract_tag,
	long long reference_events,
	long long mapped_events
)
: Statistics(run)
, match_events(0)
, used_events(0)
, oversize_events(0)
, conflict_events(0)
, detector_(detector)
, tag_(tag.empty() ? "-" : tag)
, extract_tag_(extract_tag.empty() ? "-" : extract_tag)
, reference_events_(reference_events)
, mapped_events_(mapped_events) {

}


void MatchTriggerStatistics::Write() {
	Statistics::Write<MatchTriggerStatistics>("match");
}


void MatchTriggerStatistics::Print() const {
	std::cout << detector_ << " "
		<< (tag_ == "-" ? "origin" : tag_) << " trigger"
		<< (extract_tag_.empty() ? "" : " extracting " + extract_tag_) << "\n"
		<< "Match rate "
		<< match_events << " / " << reference_events_ << "  "
		<< double(match_events) / double(reference_events_) << "\n"
		<< "Mapped events used rate "
		<< used_events << " / " << mapped_events_ << "  "
		<< double(used_events) / double(mapped_events_) << "\n"
		<< "Oversize events "
		<< oversize_events << " / " << reference_events_ << "  "
		<< double(oversize_events) / double(reference_events_) << "\n"
		<< "Conflict events "
		<< conflict_events << " / " << mapped_events_ << "  "
		<< double(conflict_events) / double(mapped_events_) << "\n";
}


std::string MatchTriggerStatistics::Title() const {
	return "run,detector,tag,extract"
		",match,used,oversize,conflict,reference,mapped"
		",match_rate,used_rate,oversize_rate,conflict_rate"
		+ title_time;
}


std::string MatchTriggerStatistics::Key() const {
	return Statistics::Key() + detector_ + tag_ + extract_tag_;
}


std::istream& operator>>(
	std::istream& is,
	MatchTriggerStatistics &statistics
) {
	CsvLineReader reader(is);
	std::string tmp;
	reader >> statistics.run_ >> statistics.detector_
		>> statistics.tag_ >> statistics.extract_tag_
		>> statistics.match_events >> statistics.used_events
		>> statistics.oversize_events >> statistics.conflict_events
		>> statistics.reference_events_ >> statistics.mapped_events_
		>> tmp >> tmp >> tmp >> tmp;
	statistics.store_time_ = reader.ReadTime();

	return is;
}


std::ostream& operator<<(
	std::ostream& os,
	const MatchTriggerStatistics &sta
) {
	os << std::setw(4) << std::setfill('0') << sta.run_
		<< "," << sta.detector_
		<< "," << sta.tag_
		<< "," << sta.extract_tag_
		<< "," << sta.match_events
		<< "," << sta.used_events
		<< "," << sta.oversize_events
		<< "," << sta.conflict_events
		<< "," << sta.reference_events_
		<< "," << sta.mapped_events_
		<< "," << double(sta.match_events) / double(sta.reference_events_)
		<< "," << double(sta.used_events) / double(sta.mapped_events_)
		<< "," << double(sta.oversize_events) / double(sta.reference_events_)
		<< "," << double(sta.conflict_events) / double(sta.mapped_events_);
	WriteStatisticsTime(os, sta.store_time_);
	return os;
}

}		// namespace ribll