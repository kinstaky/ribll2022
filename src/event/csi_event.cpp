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


void CircularCsiFundamentalEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("match", &match);
	tree->SetBranchAddress("time", time);
	tree->SetBranchAddress("energy", energy);
}


void CircularCsiFundamentalEvent::SetupOutput(TTree *tree) {
	tree->Branch("match", &match, "m/O");
	tree->Branch("time", time, "t[12]/D");
	tree->Branch("energy", energy, "e[12]/D");
}


void SquareCsiFundamentalEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("match", &match);
	tree->SetBranchAddress("time", time);
	tree->SetBranchAddress("energy", energy);
}


void SquareCsiFundamentalEvent::SetupOutput(TTree *tree) {
	tree->Branch("match", &match, "m/O");
	tree->Branch("time", time, "t[4]/D");
	tree->Branch("energy", energy, "e[4]/D");
}

}