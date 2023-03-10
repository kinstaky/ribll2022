#include "include/event/trigger_event.h"

namespace ribll {

void TriggerEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("time", &time);
	tree->SetBranchAddress("timestamp", &timestamp);
	tree->SetBranchAddress("cfd", &cfd_flag);
}


void TriggerEvent::SetupOutput(TTree *tree) {
	tree->Branch("time", &time, "t/D");
	tree->Branch("timestamp", &timestamp, "ts/L");
	tree->Branch("cfd", &cfd_flag, "cfd/O");
}

}