#include "include/statistics/track_dssd_statistics.h"

namespace ribll {

TrackDssdStatistics::TrackDssdStatistics(
	unsigned int run,
	const std::string &telescope,
	const std::string &tag
)
: Statistics(run)
, total(0)
, flag1(0)
, flag3(0)
, flag5(0)
, flag6(0)
, flag7(0)
, telescope_(telescope)
, tag_(tag.empty() ? "-" : tag) {

}


void TrackDssdStatistics::Write() {
	Statistics::Write<TrackDssdStatistics>("track-dssd");
}


void TrackDssdStatistics::Print() const {
	std::cout << "Track DSSD of " << telescope_ << " with tag " << tag_
		<< " in run " << run_ << "\n"
		<< "D1 events " << flag1 << " / " << total << "  "
		<< double(flag1) / double(total) << "\n"
		<< "D1D2 events " << flag3 << " / "  << total << "  "
		<< double(flag3) / double(total) << "\n"
		<< "D1D3 events " << flag5 << " / " << total << "  "
		<< double(flag5) / double(total) << "\n"
		<< "D2D3 events " << flag6 << " / " << total << "  "
		<< double(flag6) / double(total) << "\n"
		<< "D1D2D3 events " << flag7 << " / " << total << "  "
		<< double(flag7) / double(total) << "\n";
}


std::string TrackDssdStatistics::Title() const {
	return "run,telescope,tag"
		",flag1,flag3,flag5,flag6,flag7,total"
		",1_rate,3_rate,5_rate,6_rate,7_rate"
		+title_time;
}


std::string TrackDssdStatistics::Key() const {
	return Statistics::Key() + telescope_ + tag_;
}


std::istream& operator>>(
	std::istream &is,
	TrackDssdStatistics &statistics
) {
	CsvLineReader reader(is);
	std::string tmp;
	reader >> statistics.run_ >> statistics.telescope_ >> statistics.tag_
		>> statistics.flag1 >> statistics.flag3 >> statistics.flag5
		>> statistics.flag6 >> statistics.flag7 >> statistics.total
		>> tmp >> tmp >> tmp >> tmp >> tmp;
	statistics.store_time_ = reader.ReadTime();
	return is;
}


std::ostream& operator<<(std::ostream &os, const TrackDssdStatistics &sta) {
	os << std::setw(4) << std::setfill('0') << sta.run_
		<< "," << sta.telescope_
		<< "," << sta.tag_
		<< "," << sta.flag1
		<< "," << sta.flag3
		<< "," << sta.flag5
		<< "," << sta.flag6
		<< "," << sta.flag7
		<< "," << sta.total
		<< "," << double(sta.flag1) / double(sta.total)
		<< "," << double(sta.flag3) / double(sta.total)
		<< "," << double(sta.flag5) / double(sta.total)
		<< "," << double(sta.flag6) / double(sta.total)
		<< "," << double(sta.flag7) / double(sta.total);
	WriteStatisticsTime(os, sta.store_time_);
	return os;
}

}	// namespace ribll