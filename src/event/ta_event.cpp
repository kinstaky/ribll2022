#include "include/event/ta_event.h"

namespace ribll {

void TaEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"num").c_str(), &num);
	tree->SetBranchAddress((prefix+"flag").c_str(), flag);
	tree->SetBranchAddress((prefix+"energy").c_str(), energy);
	tree->SetBranchAddress((prefix+"radius").c_str(), radius);
	tree->SetBranchAddress((prefix+"theta").c_str(), theta);
	tree->SetBranchAddress((prefix+"phi").c_str(), phi);
}


void TaEvent::SetupOutput(TTree *tree) {
	tree->Branch("particle", &num, "num/s");
	tree->Branch("flag", flag, "f[num]/s");
	tree->Branch("energy", energy, "e[num][4]/D");
	tree->Branch("radius", radius, "r[num][4]/D");
	tree->Branch("theta", theta, "theta[num][4]/D");
	tree->Branch("phi", phi, "phi[num][4]/D");
}

}