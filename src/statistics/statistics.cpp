#include "include/statistics/statistics.h"

#include <ctime>
#include <sstream>
#include <iomanip>

namespace ribll {


Statistics::Statistics(unsigned int run)
: run_(run) {
}


std::string Statistics::Key() const {
	std::stringstream ss;
	ss << std::setw(4) << std::setfill('0') << run_;
	return ss.str();
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

}			// namespace ribll