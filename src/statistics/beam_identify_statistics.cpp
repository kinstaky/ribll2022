#include "include/statistics/beam_identify_statistics.h"

namespace ribll {

BeamIdentifyStatistics::BeamIdentifyStatistics(
	unsigned int run,
	long long total
)
: Statistics(run)
, total_(total) {
}


void BeamIdentifyStatistics::Write() {
	Statistics::Write<BeamIdentifyStatistics>("beam");
}


void BeamIdentifyStatistics::Print() const {
	std::cout << "Fitting result of 14C is "
		<< const14c << ", " << mean14c << ", " << sigma14c << "\n"
		<< "Rate of 14C is " << c14 << " / " << total_ << "  "
		<< double(c14) / double(total_) << "\n";
}


std::string BeamIdentifyStatistics::Title() const {
	return "run,type,const,mean,sigma,number,total,rate" + title_time;
}


std::string BeamIdentifyStatistics::Key() const {
	return Statistics::Key() + "14C";
}


std::istream& operator>>(
	std::istream &is,
	BeamIdentifyStatistics &statistics
) {

	CsvLineReader reader(is);
	std::string tmp;
	reader >> statistics.run_ >> tmp
		>> statistics.const14c >> statistics.mean14c >> statistics.sigma14c
		>> statistics.c14 >> statistics.total_ >> tmp;
	statistics.store_time_ = reader.ReadTime();

	return is;
}


std::ostream& operator<<(
	std::ostream &os,
	const BeamIdentifyStatistics &sta
) {
	os << std::setw(4) << std::setfill('0') << sta.run_
		<< ",14C"
		<< "," << sta.const14c
		<< "," << sta.mean14c
		<< "," << sta.sigma14c
		<< "," << sta.c14
		<< "," << sta.total_
		<< "," << double(sta.c14) / double(sta.total_);

	// write store time
	WriteStatisticsTime(os, sta.store_time_);

	return os;
}


}		// namespace ribll