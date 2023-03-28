#include "include/detector/adssd.h"

#include <iostream>


namespace ribll {

Adssd::Adssd(unsigned int run, const std::string &name)
: Dssd(run, name) {
}


int Adssd::MatchTrigger(const std::string &trigger_tag, double, double) {
	return Detector::VmeMatchTrigger<DssdFundamentalEvent>(trigger_tag);
}

}		// namespace ribll