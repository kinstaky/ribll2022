#include "include/statistics/center_statistics.h"

namespace ribll {

CenterStatistics::CenterStatistics(
	unsigned int run,
	const std::string &telescope,
	const std::string &tag
)
: Statistics(run)
, telescope_(telescope)
, tag_(tag.empty() ? "-" : tag) {

	for (size_t i = 0; i < 3; ++i) {
		x_offset[i] = 0.0;
		y_offset[i] = 0.0;
	}
}


void CenterStatistics::Write() {
	Statistics::Write<CenterStatistics>("center");
}


void CenterStatistics::Print() const {
	std::cout
		<< "x0 " << x_offset[0] << "\n"
		<< "y0 " << y_offset[0] << "\n"
		<< "x1 " << x_offset[1] << "\n"
		<< "y1 " << y_offset[1] << "\n"
		<< "x2 " << x_offset[2] << "\n"
		<< "y2 " << y_offset[2] << "\n";
}


std::string CenterStatistics::Title() const {
	return "run,telescope,tag,x0,y0,x1,y1,x2,y2" + title_time;
}


std::string CenterStatistics::Key() const {
	return Statistics::Key() + telescope_ + tag_;
}


std::istream& operator>>(
	std::istream &is,
	CenterStatistics &statistics
) {
	CsvLineReader reader(is);
	reader >> statistics.run_ >> statistics.telescope_ >> statistics.tag_
		>> statistics.x_offset[0] >> statistics.y_offset[0]
		>> statistics.x_offset[1] >> statistics.y_offset[1]
		>> statistics.x_offset[2] >> statistics.y_offset[2];
	statistics.store_time_ = reader.ReadTime();
	return is;
}


std::ostream& operator<<(
	std::ostream &os,
	const CenterStatistics &sta
) {
	os << std::setw(4) << std::setfill('0') << sta.run_
		<< "," << sta.telescope_
		<< "," << sta.tag_
		<< "," << sta.x_offset[0]
		<< "," << sta.y_offset[0]
		<< "," << sta.x_offset[1]
		<< "," << sta.y_offset[1]
		<< "," << sta.x_offset[2]
		<< "," << sta.y_offset[2];
	WriteStatisticsTime(os, sta.store_time_);
	return os;
}

}	// namespace ribll