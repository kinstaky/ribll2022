#include "include/statistics.h"

#include <ctime>
#include <sstream>
#include <iomanip>

namespace ribll {

const std::string title_time = "year,month,day,hour,minute,second,unix_time";

time_t ReadStatisticsTime(std::istream &is) {
	// temporary variable to store time
	std::string tmp;
	for (int i = 0; i < 6; ++i) std::getline(is, tmp, ',');
	// statistics time in unix time
	time_t result;
	is >> result;
	return result;
}

void WriteStatisticsTime(std::ostream &os, time_t statistics_time) {
	// output statistics time
	tm *t = localtime(&statistics_time);
	os << "," << t->tm_year+1900
		<< "," << t->tm_mon+1
		<< "," << t->tm_mday
		<< "," << t->tm_hour
		<< "," << t->tm_min
		<< "," << t->tm_sec
		<< "," << statistics_time;
}

//-----------------------------------------------------------------------------
//								Statistics
//-----------------------------------------------------------------------------

Statistics::Statistics(unsigned int run)
: run_(run) {
}


std::string Statistics::Key() const {
	std::stringstream ss;
	ss << std::setw(4) << std::setfill('0') << run_;
	return ss.str();
}

//-----------------------------------------------------------------------------
//								MapStatistics
//-----------------------------------------------------------------------------

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
	// buffer to store from getline
	std::string line;
	// stringstream to help convert to other type
	std::stringstream ss;

	// read run number
	std::getline(is, line, ',');
	ss.str(line);
	ss >> statistics.run_;

	// read crate id
	std::getline(is, line, ',');
	ss.str(line);
	ss >> statistics.crate_;

	statistics.store_time_ = ReadStatisticsTime(is);

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


//-----------------------------------------------------------------------------
//								AlignStatistics
//-----------------------------------------------------------------------------

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
	// buffer to store from getline
	std::string line;
	// stringstream to help convert to other type
	std::stringstream ss;

	// read run number
	std::getline(is, line, ',');
	ss.str(line);
	ss >> statistics.run_;

	// read align events
	std::getline(is, line, ',');
	ss.str(line);
	ss >> statistics.align_events;

	// read oversize events
	std::getline(is, line, ',');
	ss.str(line);
	ss >> statistics.oversize_events;

	// read total XIA events
	std::getline(is, line, ',');
	ss.str(line);
	ss >> statistics.xia_events_;

	// read vme events
	std::getline(is, line, ',');
	ss.str(line);
	ss >> statistics.vme_events_;

	// read XIA and VME alignment rate
	std::getline(is, line, ',');
	std::getline(is, line, ',');

	// read calibration parameters
	for (size_t i = 0; i < 2; ++i) {
		std::getline(is, line, ',');
		ss.str(line);
		ss >> statistics.calibration_param_[i];
	}

	// read store time
	statistics.store_time_ = ReadStatisticsTime(is);

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

}			// namespace ribll