#include "include/event/pd_info_event.h"

namespace ribll {

void PdInfoEvent::SetupInput(TTree *tree, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"csi_index").c_str(), &csi_index);
	tree->SetBranchAddress((prefix+"tafd_energy").c_str(), &tafd_energy);
	tree->SetBranchAddress((prefix+"taf_energy").c_str(), &taf_energy);
	tree->SetBranchAddress((prefix+"t0_energy").c_str(), &t0_energy);
	tree->SetBranchAddress((prefix+"csi_channel").c_str(), &csi_channel);
	tree->SetBranchAddress((prefix+"c_channel").c_str(), c_channel);
	tree->SetBranchAddress((prefix+"c_x").c_str(), c_x);
	tree->SetBranchAddress((prefix+"c_y").c_str(), c_y);
	tree->SetBranchAddress((prefix+"d_x").c_str(), &d_x);
	tree->SetBranchAddress((prefix+"d_y").c_str(), &d_y);
	tree->SetBranchAddress((prefix+"target_x").c_str(), &tx);
	tree->SetBranchAddress((prefix+"target_y").c_str(), &ty);
	tree->SetBranchAddress((prefix+"ppac_flag").c_str(), &ppac_flag);
	tree->SetBranchAddress((prefix+"ppac_xflag").c_str(), &ppac_xflag);
	tree->SetBranchAddress((prefix+"ppac_yflag").c_str(), &ppac_yflag);
	tree->SetBranchAddress((prefix+"ppac_x").c_str(), ppac_x);
	tree->SetBranchAddress((prefix+"ppac_y").c_str(), ppac_y);
	tree->SetBranchAddress((prefix+"ppac_x_track").c_str(), &ppac_x_track);
	tree->SetBranchAddress((prefix+"ppac_y_track").c_str(), &ppac_y_track);
	tree->SetBranchAddress((prefix+"q").c_str(), &q);
	tree->SetBranchAddress((prefix+"c15_kinetic").c_str(), &c15_kinetic);
	tree->SetBranchAddress((prefix+"c15_ex").c_str(), &c15_ex);
	tree->SetBranchAddress((prefix+"run").c_str(), &run);
	tree->SetBranchAddress((prefix+"entry").c_str(), &entry);
}


void PdInfoEvent::SetupOutput(TTree *tree) {
	tree->Branch("csi_index", &csi_index, "ci/I");
	tree->Branch("tafd_energy", &tafd_energy, "tafde/D");
	tree->Branch("taf_energy", &taf_energy, "tafe/D");
	tree->Branch("t0_energy", &t0_energy, "t0e/D");
	tree->Branch("csi_channel", &csi_channel, "csic/D");
	tree->Branch("c_channel", c_channel, "c_channel[2]/D");
	tree->Branch("c_x", c_x, "cx[2]/D");
	tree->Branch("c_y", c_y, "cy[2]/D");
	tree->Branch("d_x", &d_x, "dx/D");
	tree->Branch("d_y", &d_y, "dy/D");
	tree->Branch("target_x", &tx, "tx/D");
	tree->Branch("target_y", &ty, "ty/D");
	tree->Branch("ppac_flag", &ppac_flag, "pflag/I");
	tree->Branch("ppac_xflag", &ppac_xflag, "pxflag/s");
	tree->Branch("ppac_yflag", &ppac_yflag, "pyflag/s");
	tree->Branch("ppac_x", ppac_x, "ppacx[3]/D");
	tree->Branch("ppac_y", ppac_y, "ppacy[3]/D");
	tree->Branch("ppac_x_track", &ppac_x_track, "pxt/I");
	tree->Branch("ppac_y_track", &ppac_y_track, "pyt/I");
	tree->Branch("q", &q, "q/D");
	tree->Branch("c15_kinetic", &c15_kinetic, "ck/D");
	tree->Branch("c15_ex", &c15_ex, "cex/D");
	tree->Branch("run", &run, "run/I");
	tree->Branch("entry", &entry, "entry/L");
}

};	// ribll