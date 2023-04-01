#include "include/event/telescope_event.h"

namespace ribll {

void TelescopeEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"particle").c_str(), &particle);
	tree->SetBranchAddress((prefix+"layer").c_str(), layer);
	tree->SetBranchAddress((prefix+"flag").c_str(), flag);
	tree->SetBranchAddress((prefix+"energy").c_str(), energy);
	tree->SetBranchAddress((prefix+"radius").c_str(), radius);
	tree->SetBranchAddress((prefix+"theta").c_str(), theta);
	tree->SetBranchAddress((prefix+"phi").c_str(), phi);
}


void TelescopeEvent::SetupOutput(TTree *tree) {
	tree->Branch("particle", &particle, "p/s");
	tree->Branch("layer", layer, "l[p]/s");
	tree->Branch("flag", flag, "f[p]/s");
	tree->Branch("energy", energy, "e[p][8]/D");
	tree->Branch("radius", radius, "r[p][8]/D");
	tree->Branch("theta", theta, "theta[p][8]/D");
	tree->Branch("phi", phi, "phi[p][8]/D");
}

}