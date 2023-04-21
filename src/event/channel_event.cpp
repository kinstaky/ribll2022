#include "include/event/channel_event.h"

namespace ribll {

void ChannelEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"num").c_str(), &num);
	tree->SetBranchAddress((prefix+"charge").c_str(), charge);
	tree->SetBranchAddress((prefix+"mass").c_str(), mass);
	tree->SetBranchAddress((prefix+"energy").c_str(), energy);
	tree->SetBranchAddress((prefix+"px").c_str(), px);
	tree->SetBranchAddress((prefix+"py").c_str(), py);
	tree->SetBranchAddress((prefix+"pz").c_str(), pz);
	tree->SetBranchAddress((prefix+"r").c_str(), r);
	tree->SetBranchAddress((prefix+"theta").c_str(), theta);
	tree->SetBranchAddress((prefix+"phi").c_str(), phi);
}


void ChannelEvent::SetupOutput(TTree *tree) {
	tree->Branch("num", &num, "num/s");
	tree->Branch("charge", charge, "Z[num]/s");
	tree->Branch("mass", mass, "A[num]/s");
	tree->Branch("energy", energy, "e[num]/D");
	tree->Branch("px", px, "px[num]/D");
	tree->Branch("py", py, "py[num]/D");
	tree->Branch("pz", pz, "pz[num]/D");
	tree->Branch("r", r, "r[num]/D");
	tree->Branch("theta", theta, "theta[num]/D");
	tree->Branch("phi", phi, "phi[num]/D");
}

}	// namespace ribll