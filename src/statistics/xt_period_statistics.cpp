#include "include/statistics/xt_period_statistics.h"

namespace ribll {

XiaTriggerPeriodStatistics::XiaTriggerPeriodStatistics(
	unsigned int run,
	double period
)
: Statistics(run)
, period_(period) {
}


void XiaTriggerPeriodStatistics::Write() {
	Statistics::Write<XiaTriggerPeriodStatistics>("xt-period");
}


void XiaTriggerPeriodStatistics::Print() const {
	std::cout << "Mininum period time is " << period_ << " ns\n";
}


std::string XiaTriggerPeriodStatistics::Title() const {
	return "run,period" + title_time;
}


std::istream& operator>>(
	std::istream &is,
	XiaTriggerPeriodStatistics &statistics
) {
	CsvLineReader reader(is);
	std::string tmp;
	reader >> statistics.run_ >> statistics.period_;
	statistics.store_time_ = reader.ReadTime();

	return is;
}


std::ostream& operator<<(
	std::ostream& os,
	const XiaTriggerPeriodStatistics &sta
) {
	os << std::setw(4) << std::setfill('0') << sta.run_
		<< "," << sta.period_;

	WriteStatisticsTime(os, sta.store_time_);
	return os;
}

}		// namespace ribll