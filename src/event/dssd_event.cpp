#include "include/event/dssd_event.h"

namespace ribll {

void DssdMapEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("side", &side);
	tree->SetBranchAddress("strip", &strip);
	tree->SetBranchAddress("time", &time);
	tree->SetBranchAddress("energy", &energy);
}


void DssdMapEvent::SetupOutput(TTree *tree) {
	tree->Branch("side", &side, "side/s");
	tree->Branch("strip", &strip, "s/s");
	tree->Branch("time", &time, "t/D");
	tree->Branch("energy", &energy, "e/D");
}


void DssdFundamentalEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("front_hit", &front_hit);
	tree->SetBranchAddress("back_hit", &back_hit);
	tree->SetBranchAddress("front_strip", front_strip);
	tree->SetBranchAddress("back_strip", back_strip);
	tree->SetBranchAddress("front_time", front_time);
	tree->SetBranchAddress("back_time", back_time);
	tree->SetBranchAddress("front_energy", front_energy);
	tree->SetBranchAddress("back_energy", back_energy);
}


void DssdFundamentalEvent::SetupOutput(TTree *tree) {
	tree->Branch("front_hit", &front_hit, "fhit/s");
	tree->Branch("back_hit", &back_hit, "bhit/s");
	tree->Branch("front_strip", front_strip, "fs[fhit]/s");
	tree->Branch("back_strip", back_strip, "bs[bhit]/s");
	tree->Branch("front_time", front_time, "ft[fhit]/D");
	tree->Branch("back_time", back_time, "bt[bhit]/D");
	tree->Branch("front_energy", front_energy, "fe[fhit]/D");
	tree->Branch("back_energy", back_energy, "be[bhit]/D");
}


}	// namespace ribll