#include "include/detector/tafd.h"

#include <TMath.h>

namespace ribll {

const ROOT::Math::Polar3DVector tafd_center(135.0, 0.0, 0.0);
const std::pair<double, double> tafd_radius_ranges[6] = {
	{68, 170.5},
	{68, 170.5},
	{68, 170.5},
	{68, 170.5},
	{68, 170.5},
	{68, 170.5}
};
const std::pair<double, double> tafd_phi_ranges[6] = {
	{117.6*TMath::DegToRad(), 62.4*TMath::DegToRad()},
	{57.6*TMath::DegToRad(), 2.4*TMath::DegToRad()},
	{-2.4*TMath::DegToRad(), -57.6*TMath::DegToRad()},
	{-62.4*TMath::DegToRad(), -117.6*TMath::DegToRad()},
	{-122.4*TMath::DegToRad(), -177.6*TMath::DegToRad()},
	{177.6*TMath::DegToRad(), 122.4*TMath::DegToRad()}
};


Tafd::Tafd(unsigned int run, unsigned int index, const std::string &tag)
: Adssd(run, "tafd"+std::to_string(index), tag)
, index_(index) {

	center_ = tafd_center;
	radius_range_ = tafd_radius_ranges[index];
	phi_range_ = tafd_phi_ranges[index];
}


int Tafd::MatchTrigger(double, double) {
	if (name_ == "tafd0" || name_ == "tafd1") {
		return Detector::VmeMatchTrigger<DssdFundamentalEvent>();
	}
	std::cerr << "Error: Use ExtractTrigger instead.\n";
	return -1;
}


int Tafd::ExtractTrigger(
	double window_left,
	double window_right
) {
	if (name_ == "tafd0" || name_ == "tafd1") {
		std::cerr << "Error: Use MatchTrigger instead.\n";
		return -1;
	}
	return Dssd::ExtractTrigger(window_left, window_right);
}


int Tafd::NormalizeSides(TChain *chain, bool iteration) {
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


bool Tafd::NormEnergyCheck(size_t, const DssdFundamentalEvent &event) const {
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