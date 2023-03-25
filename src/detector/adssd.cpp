#include "include/detector/adssd.h"

#include <iostream>

#include "include/event/dssd_event.h"

namespace ribll {

Adssd::Adssd(unsigned int run, const std::string &name)
: Dssd(run, name) {
}


int Adssd::MatchTrigger(const std::string &trigger_tag, double, double) {
	return Detector::VmeMatchTrigger<DssdFundamentalEvent>(trigger_tag);
}


Taf::Taf(unsigned int run, unsigned int index)
: Adssd(run, "taf"+std::to_string(index)) {
}


int Taf::MatchTrigger(const std::string &trigger_tag, double, double) {
	if (name_ == "taf0" || name_ == "taf1") { 
		return Detector::VmeMatchTrigger<DssdFundamentalEvent>(trigger_tag);
	}
	std::cerr << "Error: Use ExtractTrigger instead.\n";
	return -1;
}


int Taf::ExtractTrigger(
	const std::string &trigger_tag,
	double window_left,
	double window_right
) {
	if (name_ == "taf0" || name_ == "taf1") {
		std::cerr << "Error: Use MatchTrigger instead.\n";
		return -1;
	}
	return Dssd::ExtractTrigger(trigger_tag, window_left, window_right);
}


Tab::Tab(unsigned int run, unsigned int index)
: Adssd(run, "tab"+std::to_string(index)) {
}


}		// namespace ribll