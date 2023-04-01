#include "include/event/tof_event.h"

namespace ribll {

void TofMapEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"time").c_str(), &time);
	tree->SetBranchAddress((prefix+"index").c_str(), &index);
	tree->SetBranchAddress((prefix+"cfd").c_str(), &cfd_flag);
}


void TofMapEvent::SetupOutput(TTree *tree) {
	tree->Branch("time", &time, "t/D");
	tree->Branch("index", &index, "index/s");
	tree->Branch("cfd", &cfd_flag, "cfd/O");
}


void TofFundamentalEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"time").c_str(), time);
	tree->SetBranchAddress((prefix+"cfd").c_str(), &cfd_flag);
}


void TofFundamentalEvent::SetupOutput(TTree *tree) {
	tree->Branch("time", time, "t[2]/D");
	tree->Branch("cfd", &cfd_flag, "cfd/s");
}

}