#include "include/statistics.h"

#include <ctime>
#include <sstream>
#include <iomanip>

namespace ribll {

const std::string title_time = ",year,month,day,hour,minute,second,unix_time";

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
	ss.clear();
	ss >> statistics.run_;

	// read crate id
	std::getline(is, line, ',');
	ss.str(line);
	ss.clear();
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
	ss.clear();
	ss >> statistics.run_;

	// read align events
	std::getline(is, line, ',');
	ss.str(line);
	ss.clear();
	ss >> statistics.align_events;

	// read oversize events
	std::getline(is, line, ',');
	ss.str(line);
	ss.clear();
	ss >> statistics.oversize_events;

	// read total XIA events
	std::getline(is, line, ',');
	ss.str(line);
	ss.clear();
	ss >> statistics.xia_events_;

	// read vme events
	std::getline(is, line, ',');
	ss.str(line);
	ss.clear();
	ss >> statistics.vme_events_;

	// read XIA and VME alignment rate
	std::getline(is, line, ',');
	std::getline(is, line, ',');

	// read calibration parameters
	for (size_t i = 0; i < 2; ++i) {
		std::getline(is, line, ',');
		ss.str(line);
		ss.clear();
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


//-----------------------------------------------------------------------------
//							MatchTriggerStatistics
//-----------------------------------------------------------------------------

MatchTriggerStatistics::MatchTriggerStatistics(
	unsigned int run,
	const std::string &detector,
	long long reference_events,
	long long mapped_events
)
: Statistics(run)
, match_events(0)
, used_events(0)
, oversize_events(0)
, conflict_events(0)
, detector_(detector)
, reference_events_(reference_events)
, mapped_events_(mapped_events) {

}


void MatchTriggerStatistics::Write() {
	Statistics::Write<MatchTriggerStatistics>("match");
}


void MatchTriggerStatistics::Print() const {
	std::cout << "Trigger match rate "
		<< match_events << " / " << reference_events_ << "  "
		<< double(match_events) / double(reference_events_) << "\n"
		<< "Mapped events used rate "
		<< used_events << " / " << mapped_events_ << "  "
		<< double(used_events) / double(mapped_events_) << "\n"
		<< "Oversize events "
		<< oversize_events << " / " << reference_events_ << "  "
		<< double(oversize_events) / double(reference_events_) << "\n"
		<< "Conflict events "
		<< conflict_events << " / " << mapped_events_ << "  "
		<< double(conflict_events) / double(mapped_events_) << "\n";
}


std::string MatchTriggerStatistics::Title() const {
	return "run,detector,match,used,oversize,conflict,reference,mapped"
		",match_rate,used_rate,oversize_rate,conflict_rate"
		+ title_time;
}


std::string MatchTriggerStatistics::Key() const {
	return Statistics::Key() + detector_;
}


std::istream& operator>>(
	std::istream& is,
	MatchTriggerStatistics &statistics
) {
	// buffer to store from getline
	std::string line;
	// stringstream to help convert to other type
	std::stringstream ss;

	// read run number
	std::getline(is, line, ',');
	ss.str(line);
	ss.clear();
	ss >> statistics.run_;

	// read detector name
	std::getline(is, line, ',');
	ss.str(line);
	ss.clear();
	ss >> statistics.detector_;

	// read match events
	std::getline(is, line, ',');
	ss.str(line);
	ss.clear();
	ss >> statistics.match_events;

	// read used events
	std::getline(is, line, ',');
	ss.str(line);
	ss.clear();
	ss >> statistics.used_events;

	// read oversize events
	std::getline(is, line, ',');
	ss.str(line);
	ss.clear();
	ss >> statistics.oversize_events;

	// read conflict events
	std::getline(is, line, ',');
	ss.str(line);
	ss.clear();
	ss >> statistics.conflict_events;

	// read reference events
	std::getline(is, line, ',');
	ss.str(line);
	ss.clear();
	ss >> statistics.reference_events_;

	// read mapped events
	std::getline(is, line, ',');
	ss.str(line);
	ss.clear();
	ss >> statistics.mapped_events_;

	// read rates
	for (int i = 0; i < 4; ++i) {
		std::getline(is, line, ',');
	}

	// read time
	statistics.store_time_ = ReadStatisticsTime(is);

	return is;
}


std::ostream& operator<<(
	std::ostream& os,
	const MatchTriggerStatistics &sta
) {
	os << std::setw(4) << std::setfill('0') << sta.run_
		<< "," << sta.detector_
		<< "," << sta.match_events
		<< "," << sta.used_events
		<< "," << sta.oversize_events
		<< "," << sta.conflict_events
		<< "," << sta.reference_events_
		<< "," << sta.mapped_events_
		<< "," << double(sta.match_events) / double(sta.reference_events_)
		<< "," << double(sta.used_events) / double(sta.mapped_events_)
		<< "," << double(sta.oversize_events) / double(sta.reference_events_)
		<< "," << double(sta.conflict_events) / double(sta.mapped_events_);

	WriteStatisticsTime(os, sta.store_time_);

	return os;
}


//-----------------------------------------------------------------------------
//							XiaTriggerPeriodStatistics
//-----------------------------------------------------------------------------

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
	// buffer to store from getline
	std::string line;
	// stringstream to help convert to other type
	std::stringstream ss;

	// read run number
	std::getline(is, line, ',');
	ss.str(line);
	ss.clear();
	ss >> statistics.run_;

	// read period
	std::getline(is, line, ',');
	ss.str(line);
	ss.clear();
	ss >> statistics.period_;

		// read time
	statistics.store_time_ = ReadStatisticsTime(is);

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



}			// namespace ribll