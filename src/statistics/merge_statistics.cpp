#include "include/statistics/merge_statistics.h"

namespace ribll {

MergeStatistics::MergeStatistics(
	unsigned int run,
	const std::string &detector,
	const std::string &tag
)
: Statistics(run)
, total(0)
, merged(0)
, one_hit(0)
, two_hit(0)
, three_hit(0)
, detector_(detector)
, tag_(tag.empty() ? "-" : tag) {
}


void MergeStatistics::Write() {
	Statistics::Write<MergeStatistics>("merge");
}


void MergeStatistics::Print() const {
	std::cout << detector_ << " "
		<< (tag_ == "-" ? "origin" : tag_) << " trigger\n"
		<< "Merge rate "
		<< merged << " / " << total << "  "
		<< double(merged) / double(total) << "\n"
		<< "1 particle rate "
		<< one_hit << " / " << total << "  "
		<< double(one_hit) / double(total) << "\n"
		<< "2 particle rate "
		<< two_hit << " / " << total << "  "
		<< double(two_hit) / double(total) << "\n"
		<< "3 particle rate "
		<< three_hit << " / " << total << "  "
		<< double(three_hit) / double(total) << "\n";
}


std::string MergeStatistics::Title() const {
	return "run,detector,tag,total,merge,merge_rate"
		",one,two,three,one_rate,two_rate,three_rate"
		+ title_time;
}


std::string MergeStatistics::Key() const {
	return Statistics::Key() + detector_ + tag_;
}


std::istream& operator>>(
	std::istream &is,
	MergeStatistics &statistics
) {
	CsvLineReader reader(is);
	std::string tmp;
	reader >> statistics.run_ >> statistics.detector_ >> statistics.tag_
		>> statistics.total >> statistics.merged >> tmp
		>> statistics.one_hit >> statistics.two_hit >> statistics.three_hit
		>> tmp >> tmp >> tmp;
	statistics.store_time_ = reader.ReadTime();
	return is;
}


std::ostream& operator<<(
	std::ostream &os,
	const MergeStatistics &sta
) {
	os << std::setw(4) << std::setfill('0') << sta.run_
		<< "," << sta.detector_
		<< "," << sta.tag_
		<< "," << sta.total
		<< "," << sta.merged
		<< "," << double(sta.merged) / double(sta.total)
		<< "," << sta.one_hit
		<< "," << sta.two_hit
		<< "," << sta.three_hit
		<< "," << double(sta.one_hit) / double(sta.total)
		<< "," << double(sta.two_hit) / double(sta.total)
		<< "," << double(sta.three_hit) / double(sta.total);
	WriteStatisticsTime(os, sta.store_time_);
	return os;
}

}	// namespace ribll