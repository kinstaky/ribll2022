#include "include/statistics/track_statistics.h"

namespace ribll {

TrackStatistics::TrackStatistics(
	unsigned int run,
	const std::string &telescope,
	const std::string &tag
)
: Statistics(run)
, total(0)
, particle1(0)
, particle2(0)
, particle3(0)
, telescope_(telescope)
, tag_(tag.empty() ? "-" : tag) {
}


void TrackStatistics::Write() {
	Statistics::Write<TrackStatistics>("track");
}


void TrackStatistics::Print() const {
	std::cout << "Track " << telescope_ << " with tag " << tag_
		<< " in run " << run_ << "\n"
		<< "1 particle events " << particle1 << " / " << total << "  "
		<< double(particle1) / double(total) << "\n"
		<< "2 particle events " << particle2 << " / " << total << "  "
		<< double(particle2) / double(total) << "\n"
		<< "3 particle events " << particle3 << " / " << total << "  "
		<< double(particle3) / double(total) << "\n";
}


std::string TrackStatistics::Title() const {
	return "run,telescope,tag,1_particle,2_particle,3_particle,total"
		",1p_rate,2p_rate,3p_rate" + title_time;
}


std::string TrackStatistics::Key() const {
	return Statistics::Key() + telescope_ + tag_;
}


std::istream& operator>>(std::istream &is, TrackStatistics &statistics) {
	CsvLineReader reader(is);
	std::string tmp;
	reader >> statistics.run_ >> statistics.telescope_ >> statistics.tag_
		>> statistics.particle1 >> statistics.particle2
		>> statistics.particle3 >> statistics.total
		>> tmp >> tmp >> tmp;
	statistics.store_time_ = reader.ReadTime();
	return is;
}


std::ostream& operator<<(
	std::ostream &os,
	const TrackStatistics &sta
) {
	os << std::setw(4) << std::setfill('0') << sta.run_
		<< "," << sta.telescope_
		<< "," << sta.tag_
		<< "," << sta.particle1
		<< "," << sta.particle2
		<< "," << sta.particle3
		<< "," << sta.total
		<< "," << double(sta.particle1) / double(sta.total)
		<< "," << double(sta.particle2) / double(sta.total)
		<< "," << double(sta.particle3) / double(sta.total);
	WriteStatisticsTime(os, sta.store_time_);
	return os;
}

}		// namespace ribll