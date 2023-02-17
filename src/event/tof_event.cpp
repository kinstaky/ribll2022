#include "include/event/tof_event.h"

namespace ribll {

void TofMapEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("time", &time);
	tree->SetBranchAddress("timestamp", &timestamp);
	tree->SetBranchAddress("index", &index);
}


void TofMapEvent::SetupOutput(TTree *tree) {
	tree->Branch("time", &time, "t/D");
	tree->Branch("timestamp", &timestamp, "ts/L");
	tree->Branch("index", &index, "index/s");
}


void TofFundamentalEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("time", time);
}


void TofFundamentalEvent::SetupOutput(TTree *tree) {
	tree->Branch("time", time, "t[2]/D");
}

}