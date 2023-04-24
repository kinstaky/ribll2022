#include "include/statistics/normalize_statistics.h"

namespace ribll {

NormalizeStatistics::NormalizeStatistics(
	unsigned int run,
	const std::string &detector,
	const std::string &tag,
	unsigned int end_run,
	int iteration
)
: Statistics(run)
, detector_(detector)
, tag_(tag.empty() ? "-" : tag)
, end_run_(end_run)
, iteration_(iteration) {
}


void NormalizeStatistics::Write() {
	Statistics::Write<NormalizeStatistics>("normalize");
}


void NormalizeStatistics::Print() const {
	std::cout << detector_ << " "
		<< (tag_ == "-" ? "origin" : tag_) << " trigger\n"
		<< "Normalize " << run_ << " to " << end_run_
		<< " in iteration " << iteration_ << "\n";
}


std::string NormalizeStatistics::Title() const {
	return "run,detector,tag,end_run,iteration"+ title_time;
}


std::string NormalizeStatistics::Key() const {
	return Statistics::Key() + detector_ + tag_ + std::to_string(iteration_);
}


std::istream& operator>>(
	std::istream &is,
	NormalizeStatistics &statistics
) {
	CsvLineReader reader(is);
	std::string tmp;
	reader >> statistics.run_ >> statistics.detector_ >> statistics.tag_
		>> statistics.end_run_ >> statistics.iteration_;
	statistics.store_time_ = reader.ReadTime();
	return is;
}


std::ostream& operator<<(
	std::ostream &os,
	const NormalizeStatistics &sta
) {
	os << std::setw(4) << std::setfill('0') << sta.run_
		<< "," << sta.detector_
		<< "," << sta.tag_
		<< "," << sta.end_run_
		<< "," << sta.iteration_;
	WriteStatisticsTime(os, sta.store_time_);
	return os;
}

}	// namespace ribll