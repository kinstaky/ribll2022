#include "include/event/ssd_event.h"

namespace ribll {

void SsdEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("cfd", &cfd_flag);
	tree->SetBranchAddress("time", &time);
	tree->SetBranchAddress("energy", &energy);
}


void SsdEvent::SetupOutput(TTree *tree) {
	tree->Branch("cfd", &cfd_flag, "cfd/O");
	tree->Branch("time", &time, "t/D");
	tree->Branch("energy", &energy, "e/D");
}

}		// namespace ribll