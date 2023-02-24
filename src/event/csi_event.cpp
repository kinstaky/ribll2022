#include "include/event/csi_event.h"


namespace ribll {

void CsiMapEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("index", &index);
	tree->SetBranchAddress("time", &time);
	tree->SetBranchAddress("energy", &energy);
}


void CsiMapEvent::SetupOutput(TTree *tree) {
	tree->Branch("index", &index, "index/s");
	tree->Branch("time", &time, "t/D");
	tree->Branch("energy", &energy, "e/D");
}


void CsiFundamentalEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("time", &time);
	tree->SetBranchAddress("energy", &energy);
}


void CsiFundamentalEvent::SetupOutput(TTree *tree) {
	tree->Branch("time", &time, "t/D");
	tree->Branch("energy", &energy, "e/D");
}

};