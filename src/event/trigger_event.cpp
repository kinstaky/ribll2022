#include "include/event/trigger_event.h"

namespace ribll {

void TriggerEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("time", &time);
}


void TriggerEvent::SetupOutput(TTree *tree) {
	tree->Branch("time", &time, "t/D");
}

}