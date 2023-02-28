#include "include/detector/adssd.h"

#include <iostream>

namespace ribll {

Adssd::Adssd(unsigned int run, const std::string &name)
: Dssd(run, name) {
}


Taf::Taf(unsigned int run, unsigned int index)
: Adssd(run, "taf"+std::to_string(index)) {
}


int Taf::MatchTrigger(double window_left, double window_right) {
	if (name_ == "taf0" || name_ == "taf1") {
		std::cerr << "Error: MatchTrigger with vt instead of "
			<< name_ << "\n";
		return -1;
	}
	return Dssd::MatchTrigger(window_left, window_right);
}

Tab::Tab(unsigned int run, unsigned int index)
: Adssd(run, "tab"+std::to_string(index)) {
}


int Tab::MatchTrigger(double, double) {
	std::cerr << "Error: MatchTrigger with vt instead of "
		<< name_ << "\n";
	return -1;
}

}