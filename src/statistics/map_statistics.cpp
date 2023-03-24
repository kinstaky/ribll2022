#include "include/statistics/map_statistics.h"

namespace ribll {

MapStatistics::MapStatistics(unsigned int run, unsigned int crate)
: Statistics(run)
, crate_(crate) {
}


void MapStatistics::Write() {
	Statistics::Write<MapStatistics>("map");
}


void MapStatistics::Print() const {
}


std::string MapStatistics::Title() const {
	return "run,crate" + title_time;
}


std::string MapStatistics::Key() const {
	return Statistics::Key() + std::to_string(crate_);
}


std::istream& operator>>(
	std::istream &is,
	MapStatistics &statistics
) {
	CsvLineReader reader(is);
	std::string tmp;
	reader >> statistics.run_ >> statistics.crate_;
	statistics.store_time_ = reader.ReadTime();

	return is;
}


std::ostream& operator<<(
	std::ostream &os,
	const MapStatistics &statistics
) {
	os << std::setw(4) << std::setfill('0') << statistics.run_
		<< "," << statistics.crate_;

	WriteStatisticsTime(os, statistics.store_time_);

	return os;
}

}		// namespace ribll