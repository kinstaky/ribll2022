#include "include/statistics/align_statistics.h"

namespace ribll {

AlignStatistics::AlignStatistics(
	unsigned int run,
	long long xia_events,
	long long vme_events,
	double *calibration_parameters
)
: Statistics(run)
, align_events(0)
, oversize_events(0)
, xia_events_(xia_events)
, vme_events_(vme_events) {

	calibration_param_[0] = calibration_parameters[0];
	calibration_param_[1] = calibration_parameters[1];
}


void AlignStatistics::Write() {
	Statistics::Write<AlignStatistics>("align");
}


void AlignStatistics::Print() const {
	std::cout << "xia alignment rate "
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
	return "run,align,oversize,xia,vme,xia_rate,vme_rate,p0,p1" + title_time;
}


std::istream& operator>>(
	std::istream &is,
	AlignStatistics &statistics
) {
	CsvLineReader reader(is);
	std::string tmp;
	reader >> statistics.run_
		>> statistics.align_events >> statistics.oversize_events
		>> statistics.xia_events_ >> statistics.vme_events_
		>> tmp >> tmp
		>> statistics.calibration_param_[0]
		>> statistics.calibration_param_[1];
	statistics.store_time_ = reader.ReadTime();

	return is;
}


std::ostream& operator<<(
	std::ostream &os,
	const AlignStatistics &sta
) {
	os << std::setw(4) << std::setfill('0') << sta.run_
		<< "," << sta.align_events
		<< "," << sta.oversize_events
		<< "," << sta.xia_events_
		<< "," << sta.vme_events_
		<< "," << double(sta.align_events) / double(sta.xia_events_)
		<< "," << double(sta.align_events) / double(sta.vme_events_)
		<< "," << std::setprecision(15) << sta.calibration_param_[0]
		<< "," << std::setprecision(15) << sta.calibration_param_[1];
	WriteStatisticsTime(os, sta.store_time_);
	return os;
}

}		// namespace ribll