#include "include/event/csi_event.h"


namespace ribll {

void CsiMapEvent::SetupInput(TTree *tree, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"index").c_str(), &index);
	tree->SetBranchAddress((prefix+"cfd").c_str(), &cfd_flag);
	tree->SetBranchAddress((prefix+"time").c_str(), &time);
	tree->SetBranchAddress((prefix+"energy").c_str(), &energy);
}


void CsiMapEvent::SetupOutput(TTree *tree) {
	tree->Branch("index", &index, "index/s");
	tree->Branch("cfd", &cfd_flag, "cfd/O");
	tree->Branch("time", &time, "t/D");
	tree->Branch("energy", &energy, "e/D");
}


void CircularCsiFundamentalEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"match").c_str(), &match);
	tree->SetBranchAddress((prefix+"cfd").c_str(), &cfd_flag);
	tree->SetBranchAddress((prefix+"time").c_str(), time);
	tree->SetBranchAddress((prefix+"energy").c_str(), energy);
}


void CircularCsiFundamentalEvent::SetupOutput(TTree *tree) {
	tree->Branch("match", &match, "m/O");
	tree->Branch("cfd", &cfd_flag, "cfd/s");
	tree->Branch("time", time, "t[12]/D");
	tree->Branch("energy", energy, "e[12]/D");
}


void SquareCsiFundamentalEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"match").c_str(), &match);
	tree->SetBranchAddress((prefix+"cfd").c_str(), &cfd_flag);
	tree->SetBranchAddress((prefix+"time").c_str(), time);
	tree->SetBranchAddress((prefix+"energy").c_str(), energy);
}


void SquareCsiFundamentalEvent::SetupOutput(TTree *tree) {
	tree->Branch("match", &match, "m/O");
	tree->Branch("cfd", &cfd_flag, "cfd/s");
	tree->Branch("time", time, "t[4]/D");
	tree->Branch("energy", energy, "e[4]/D");
}

}