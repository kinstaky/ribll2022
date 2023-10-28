#include "include/event/ta_event.h"

namespace ribll {

void TaEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"num").c_str(), &num);
	tree->SetBranchAddress((prefix+"flag").c_str(), flag);
	tree->SetBranchAddress((prefix+"energy").c_str(), energy);
	tree->SetBranchAddress((prefix+"time").c_str(), time);
	tree->SetBranchAddress((prefix+"radius").c_str(), radius);
	tree->SetBranchAddress((prefix+"theta").c_str(), theta);
	tree->SetBranchAddress((prefix+"phi").c_str(), phi);
	tree->SetBranchAddress((prefix+"front_strip").c_str(), front_strip);
	tree->SetBranchAddress((prefix+"back_strip").c_str(), back_strip);
}


void TaEvent::SetupOutput(TTree *tree) {
	tree->Branch("num", &num, "num/I");
	tree->Branch("flag", flag, "f[num]/s");
	tree->Branch("energy", energy, "e[num][4]/D");
	tree->Branch("time", time, "t[num][4]/D");
	tree->Branch("radius", radius, "r[num]/D");
	tree->Branch("theta", theta, "theta[num]/D");
	tree->Branch("phi", phi, "phi[num]/D");
	tree->Branch("front_strip", front_strip, "fs[num]/s");
	tree->Branch("back_strip", back_strip, "bs[num]/s");
}

}