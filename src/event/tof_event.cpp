#include "include/event/tof_event.h"

namespace ribll {

void TofMapEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("time", &time);
	tree->SetBranchAddress("index", &index);
	tree->SetBranchAddress("cfd", &cfd_flag);
}


void TofMapEvent::SetupOutput(TTree *tree) {
	tree->Branch("time", &time, "t/D");
	tree->Branch("index", &index, "index/s");
	tree->Branch("cfd", &cfd_flag, "cfd/O");
}


void TofFundamentalEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("time", time);
	tree->SetBranchAddress("cfd", &cfd_flag);
}


void TofFundamentalEvent::SetupOutput(TTree *tree) {
	tree->Branch("time", time, "t[2]/D");
	tree->Branch("cfd", &cfd_flag, "cfd/s");
}

}