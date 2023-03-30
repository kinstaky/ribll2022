#include "include/statistics/t0_hit_statistics.h"

namespace ribll {

T0HitStatistics::T0HitStatistics(
	unsigned int run,
	const std::string &tag,
	long long total
)
: Statistics(run)
, d1_single_hit(0)
, d1_multi_hit(0)
, d2_single_hit(0)
, d2_multi_hit(0)
, d1d2_single_hit(0)
, d1d2_multi_hit(0)
, tag_(tag.empty() ? "-" : tag)
, total_(total) {
}

void T0HitStatistics::Write() {
	Statistics::Write<T0HitStatistics>("t0hit");
}


void T0HitStatistics::Print() const {
	std::cout << "Hit information of " << tag_ << " t0 in " << run_ << "\n"
		<< "t0d1 hit " << d1_single_hit << " / " << d1_multi_hit
		<< "  " << double(d1_single_hit) / double(d1_multi_hit) << "\n"
		<< "t0d2 hit " << d2_single_hit << " / " << d2_multi_hit
		<< "  " << double(d2_single_hit) / double(d2_multi_hit) << "\n"
		<< "t0d1d2 hit " << d1d2_single_hit << " / " << d1d2_multi_hit
		<< "  " << double(d1d2_single_hit) / double(d1d2_multi_hit) << "\n";
}


std::string T0HitStatistics::Title() const {
	return "run,tag"
		",d1s,d1m,d2s,d2m,d1d2s,d1d2m,total"
		",d1sm_rate,d2sm_rate,d1d2sm_rate"
		+ title_time;
}


std::string T0HitStatistics::Key() const {
	return Statistics::Key() + tag_;
}


std::istream& operator>>(std::istream &is, T0HitStatistics &statistics) {
	CsvLineReader reader(is);
	std::string tmp;
	reader >> statistics.run_ >> statistics.tag_
		>> statistics.d1_single_hit >> statistics.d1_multi_hit
		>> statistics.d2_single_hit >> statistics.d2_multi_hit
		>> statistics.d1d2_single_hit >> statistics.d1d2_multi_hit
		>> statistics.total_
		>> tmp >> tmp >> tmp;
	statistics.store_time_ = reader.ReadTime();
	return is;
}


std::ostream& operator<<(std::ostream &os, const T0HitStatistics &sta) {
	os << std::setw(4) << std::setfill('0') << sta.run_
		<< "," << sta.tag_
		<< "," << sta.d1_single_hit
		<< "," << sta.d1_multi_hit
		<< "," << sta.d2_single_hit
		<< "," << sta.d2_multi_hit
		<< "," << sta.d1d2_single_hit
		<< "," << sta.d1d2_multi_hit
		<< "," << sta.total_
		<< "," << double(sta.d1_single_hit) / double(sta.d1_multi_hit)
		<< "," << double(sta.d2_single_hit) / double(sta.d2_multi_hit)
		<< "," << double(sta.d1d2_single_hit) / double(sta.d1d2_multi_hit);

	WriteStatisticsTime(os, sta.store_time_);
	return os;
}

}		// namespace ribll