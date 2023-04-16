#include "include/event/particle_type_event.h"

namespace ribll {

void ParticleTypeEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"num").c_str(), &num);
	tree->SetBranchAddress((prefix+"charge").c_str(), charge);
	tree->SetBranchAddress((prefix+"mass").c_str(), mass);
	tree->SetBranchAddress((prefix+"layer").c_str(), layer);
}


void ParticleTypeEvent::SetupOutput(TTree *tree) {
	tree->Branch("num", &num, "num/s");
	tree->Branch("charge", charge, "z[num]/s");
	tree->Branch("mass", mass, "a[num]/s");
	tree->Branch("layer", layer, "l[num]/S");
}

}	// namespace ribll