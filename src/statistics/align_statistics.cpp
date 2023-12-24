#include "include/statistics/align_statistics.h"

namespace ribll {

AlignStatistics::AlignStatistics(
	unsigned int run,
	long long xia_events,
	long long vme_events,
	unsigned int group
)
: Statistics(run)
, first_align_events(0)
, second_align_events(0)
, third_align_events(0)
, align_events(0)
, oversize_events(0)
, xia_events_(xia_events)
, vme_events_(vme_events)
, group_(group) {
}


void AlignStatistics::Write() {
	Statistics::Write<AlignStatistics>("align");
}


void AlignStatistics::Print() const {
	std::cout
		<< "first alignment rate "
		<< first_align_events << " / " << xia_events_ << " "
		<< double(first_align_events) / double(xia_events_) << "\n"
		<< "second alignment rate "
		<< second_align_events << " / " << xia_events_ << " "
		<< double(second_align_events) / double(xia_events_) << "\n"
		<< "third alignment rate "
		<< third_align_events << " / " << xia_events_ << " "
		<< double(third_align_events) / double(xia_events_) << "\n"
		<< "xia alignment rate "
		<< align_events << " / " << xia_events_ << " "
		<< double(align_events) / double(xia_events_) << "\n"
		<< "vme alignment rate "
		<< align_events << " / " << vme_events_ << " "
		<< double(align_events) / double(vme_events_) << "\n"
		<< "oversize events "
		<< oversize_events << " / " << vme_events_ << " "
		<< double(oversize_events) / double(vme_events_) << "\n";
}


std::string AlignStatistics::Title() const {
	return "run,group,align,oversize,xia,vme,xia_rate,vme_rate" + title_time;
}


std::istream& operator>>(
	std::istream &is,
	AlignStatistics &statistics
) {
	CsvLineReader reader(is);
	std::string tmp;
	reader >> statistics.run_ >> statistics.group_
		>> statistics.align_events >> statistics.oversize_events
		>> statistics.xia_events_ >> statistics.vme_events_
		>> tmp >> tmp;
	statistics.store_time_ = reader.ReadTime();

	return is;
}


std::ostream& operator<<(
	std::ostream &os,
	const AlignStatistics &sta
) {
	os << std::setw(4) << std::setfill('0') << sta.run_
		<< "," << sta.group_
		<< "," << sta.align_events
		<< "," << sta.oversize_events
		<< "," << sta.xia_events_
		<< "," << sta.vme_events_
		<< "," << double(sta.align_events) / double(sta.xia_events_)
		<< "," << double(sta.align_events) / double(sta.vme_events_);
	WriteStatisticsTime(os, sta.store_time_);
	return os;
}

}		// namespace ribll