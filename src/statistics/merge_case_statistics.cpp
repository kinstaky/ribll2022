#include "include/statistics/merge_case_statistics.h"

namespace ribll {

MergeCaseStatistics::MergeCaseStatistics(
	unsigned int run,
	const std::string &detector,
	const std::string &tag,
	const std::string &case_name,
	unsigned int minor_case,
	double tolerance
)
: Statistics(run)
, total(0)
, merged(0)
, detector_(detector)
, tag_(tag.empty() ? "-" : tag)
, case_name_(case_name)
, minor_case_(minor_case)
, tolerance_(tolerance) {
}


void MergeCaseStatistics::Write() {
	Statistics::Write<MergeCaseStatistics>("merge-case");
}


void MergeCaseStatistics::Print() const {
	std::cout << detector_ << " "
		<< (tag_ == "-" ? "origin" : tag_) << " trigger\n"
		<< "Merge in case " << case_name_ << " " << minor_case_
		<< " tolerance " << tolerance_ << "\n"
		<< "Merge rate " << merged  << " / " << total << "  "
		<< double(merged) / double(total) << "\n";
}


std::string MergeCaseStatistics::Title() const {
	return "run,detector,tag,case,minor,tolerance,total,merge,rate"+title_time;
}


std::string MergeCaseStatistics::Key() const {
	return Statistics::Key() + detector_ + tag_
		+ case_name_ + std::to_string(minor_case_);
}


std::istream& operator>>(
	std::istream &is,
	MergeCaseStatistics &statistics
) {
	CsvLineReader reader(is);
	std::string tmp;
	reader >> statistics.run_ >> statistics.detector_ >> statistics.tag_
		>> statistics.case_name_ >> statistics.minor_case_
		>> statistics.tolerance_
		>> statistics.total >> statistics.merged >> tmp;
	statistics.store_time_ = reader.ReadTime();
	return is;
}


std::ostream& operator<<(std::ostream &os, const MergeCaseStatistics &sta) {
	os << std::setw(4) << std::setfill('0') << sta.run_
		<< "," << sta.detector_
		<< "," << sta.tag_
		<< "," << sta.case_name_
		<< "," << sta.minor_case_
		<< "," << sta.tolerance_
		<< "," << sta.total
		<< "," << sta.merged
		<< "," << double(sta.merged) / double(sta.total);
	WriteStatisticsTime(os, sta.store_time_);
	return os;
}


}	// namespace ribll