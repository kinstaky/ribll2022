#include "include/detector/taf.h"

namespace ribll {

Taf::Taf(unsigned int run, unsigned int index)
: Adssd(run, "taf"+std::to_string(index))
, index_(index) {
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


int Taf::NormalizeSides(TChain *chain, bool iteration) {
	if (SideNormalize(chain, 0, 4, iteration)) {
		std::cerr << "Error: Normalize first side failed.\n";
		return -1;
	}
	if (SideNormalize(chain, 1, 1, iteration)) {
		std::cerr << "Error: Normalize second side failed.\n";
		return -1;
	}
	return 0;
}


bool Taf::NormEnergyCheck(size_t, const DssdFundamentalEvent &event) {
	if (index_ == 0) {
		if (event.front_energy[0] > 1e4 || event.back_energy[0] > 1e4) {
			return false;
		}
	} else if (index_ == 1) {
		if (event.front_energy[0] > 1e4 || event.back_energy[0] > 1e4) {
			return false;
		}
	} else {
		return false;
	}
	return true;
}

}