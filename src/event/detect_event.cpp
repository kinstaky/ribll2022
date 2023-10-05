#include "include/event/detect_event.h"

namespace ribll {

void DetectEvent::SetupInput(TTree *tree, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"taf_layer").c_str(), &taf_layer);
	tree->SetBranchAddress(
		(prefix+"taf_lost_energy").c_str(), taf_lost_energy
	);
	tree->SetBranchAddress((prefix+"taf_energy").c_str(), taf_energy);
	tree->SetBranchAddress((prefix+"t0_layer").c_str(), t0_layer);
	tree->SetBranchAddress((prefix+"t0_lost_energy").c_str(), t0_lost_energy);
	tree->SetBranchAddress((prefix+"t0_energy").c_str(), t0_energy);
	tree->SetBranchAddress((prefix+"tafx").c_str(), &tafx);
	tree->SetBranchAddress((prefix+"tafy").c_str(), &tafy);
	tree->SetBranchAddress((prefix+"tafz").c_str(), &tafz);
	tree->SetBranchAddress((prefix+"tafr").c_str(), &tafr);
	tree->SetBranchAddress((prefix+"t0x").c_str(), t0x);
	tree->SetBranchAddress((prefix+"t0y").c_str(), t0y);
	tree->SetBranchAddress((prefix+"t0z").c_str(), t0z);
	tree->SetBranchAddress((prefix+"t0r").c_str(), t0r);
	tree->SetBranchAddress((prefix+"ppacx").c_str(), ppacx);
	tree->SetBranchAddress((prefix+"ppacy").c_str(), ppacy);
	tree->SetBranchAddress((prefix+"tx").c_str(), &tx);
	tree->SetBranchAddress((prefix+"ty").c_str(), &ty);
	tree->SetBranchAddress((prefix+"valid").c_str(), &valid);
	tree->SetBranchAddress((prefix+"q").c_str(), &q);
}


void DetectEvent::SetupOutput(TTree *tree) {
	tree->Branch("taf_layer", &taf_layer, "taflayer/I");
	tree->Branch("taf_lost_energy", taf_lost_energy, "tafle[2]/D");
	tree->Branch("taf_energy", taf_energy, "tafe[2]/D");
	tree->Branch("t0_layer", t0_layer, "t0layer[2]/I");
	tree->Branch("t0_lost_energy", t0_lost_energy, "t0le[2][7]/D");
	tree->Branch("t0_energy", t0_energy, "t0e[2][7]/D");
	tree->Branch("tafx", &tafx, "tafx/D");
	tree->Branch("tafy", &tafy, "tafy/D");
	tree->Branch("tafz", &tafz, "tafz/D");
	tree->Branch("tafr", &tafr, "tafr/D");
	tree->Branch("t0x", t0x, "t0x[2][3]/D");
	tree->Branch("t0y", t0y, "t0y[2][3]/D");
	tree->Branch("t0z", t0z, "t0z[2][3]/D");
	tree->Branch("t0r", t0r, "t0r[2][3]/D");
	tree->Branch("ppacx", ppacx, "ppacx[3]/D");
	tree->Branch("ppacy", ppacy, "ppacy[3]/D");
	tree->Branch("tx", &tx, "tx/D");
	tree->Branch("ty", &ty, "ty/D");
	tree->Branch("valid", &valid, "valid/I");
	tree->Branch("q", &q, "q/D");
}


}	// ribll