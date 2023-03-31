#include "include/detector/tabd.h"

namespace ribll {

const ROOT::Math::Polar3DVector tabd_center(0.0, ROOT::Math::Pi(), 0.0);
const std::pair<double, double> tabd_radius_ranges[6] = {
	{0.068, 0.1705},
	{0.068, 0.1705},
	{0.068, 0.1705},
	{0.068, 0.1705},
	{0.068, 0.1705},
	{0.068, 0.1705}
};
const std::pair<double, double> tabd_phi_ranges[6] = {
	{62.4*TMath::DegToRad(), 117.6*TMath::DegToRad()},
	{122.4*TMath::DegToRad(), 177.6*TMath::DegToRad()},
	{-177.6*TMath::DegToRad(), -122.4*TMath::DegToRad()},
	{-117.6*TMath::DegToRad(), -62.4*TMath::DegToRad()},
	{-57.6*TMath::DegToRad(), -2.4*TMath::DegToRad()},
	{2.4*TMath::DegToRad(), 57.6*TMath::DegToRad()}
};

Tabd::Tabd(unsigned int run, unsigned int index, const std::string &tag)
: Adssd(run, "tab"+std::to_string(index), tag) {
}


int Tabd::MatchTrigger(double, double) {
	return Detector::VmeMatchTrigger<DssdFundamentalEvent>();
}


}		// namespace ribll