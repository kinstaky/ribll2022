#include "include/event/trigger_event.h"

namespace ribll {

void TriggerEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"time").c_str(), &time);
	tree->SetBranchAddress((prefix+"timestamp").c_str(), &timestamp);
	tree->SetBranchAddress((prefix+"cfd").c_str(), &cfd_flag);
}


void TriggerEvent::SetupOutput(TTree *tree) {
	tree->Branch("time", &time, "t/D");
	tree->Branch("timestamp", &timestamp, "ts/L");
	tree->Branch("cfd", &cfd_flag, "cfd/O");
}

}