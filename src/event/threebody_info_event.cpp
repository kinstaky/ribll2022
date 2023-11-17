#include "include/event/threebody_info_event.h"

namespace ribll {

void ThreeBodyInfoEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"csi_index").c_str(), &csi_index);
	tree->SetBranchAddress((prefix+"layer").c_str(), layer);
	tree->SetBranchAddress((prefix+"tafd_energy").c_str(), &tafd_energy);
	tree->SetBranchAddress((prefix+"t0_energy").c_str(), t0_energy);
	tree->SetBranchAddress((prefix+"csi_channel").c_str(), &csi_channel);
	tree->SetBranchAddress((prefix+"d1_channel").c_str(), d1_channel);
	tree->SetBranchAddress((prefix+"d2_channel").c_str(), d2_channel);
	tree->SetBranchAddress((prefix+"d3_channel").c_str(), d3_channel);
	tree->SetBranchAddress((prefix+"ssd_channel").c_str(), ssd_channel);
	tree->SetBranchAddress((prefix+"d1x").c_str(), d1x);
	tree->SetBranchAddress((prefix+"d1y").c_str(), d1y);
	tree->SetBranchAddress((prefix+"d2x").c_str(), d2x);
	tree->SetBranchAddress((prefix+"d2y").c_str(), d2y);
	tree->SetBranchAddress((prefix+"d3x").c_str(), d3x);
	tree->SetBranchAddress((prefix+"d3y").c_str(), d3y);
	tree->SetBranchAddress((prefix+"rx").c_str(), &rx);
	tree->SetBranchAddress((prefix+"ry").c_str(), &ry);
	tree->SetBranchAddress((prefix+"tx").c_str(), &tx);
	tree->SetBranchAddress((prefix+"ty").c_str(), &ty);
	tree->SetBranchAddress((prefix+"ppac_xflag").c_str(), &ppac_xflag);
	tree->SetBranchAddress((prefix+"ppac_yflag").c_str(), &ppac_yflag);
	tree->SetBranchAddress((prefix+"ppac_x").c_str(), ppac_x);
	tree->SetBranchAddress((prefix+"ppac_y").c_str(), ppac_y);
	tree->SetBranchAddress((prefix+"run").c_str(), &run);
}


void ThreeBodyInfoEvent::SetupOutput(TTree *tree) {
	tree->Branch("csi_index", &csi_index, "ci/I");
	tree->Branch("layer", layer, "layer[2]/I");
	tree->Branch("tafd_energy", &tafd_energy, "tafde/D");
	tree->Branch("t0_energy", t0_energy, "t0e[2]/D");
	tree->Branch("csi_channel", &csi_channel, "csic/D");
	tree->Branch("d1_channel", d1_channel, "d1c[2]/D");
	tree->Branch("d2_channel", d2_channel, "d2c[2]/D");
	tree->Branch("d3_channel", d3_channel, "d3c[2]/D");
	tree->Branch("ssd_channel", ssd_channel, "ssdc[3]/D");
	tree->Branch("d1x", d1x, "d1x[2]/D");
	tree->Branch("d1y", d1y, "d1y[2]/D");
	tree->Branch("d2x", d2x, "d2x[2]/D");
	tree->Branch("d2y", d2y, "d2y[2]/D");
	tree->Branch("d3x", d3x, "d3x[2]/D");
	tree->Branch("d3y", d3y, "d3y[2]/D");
	tree->Branch("rx", &rx, "rx/D");
	tree->Branch("ry", &ry, "ry/D");
	tree->Branch("tx", &tx, "tx/D");
	tree->Branch("ty", &ty, "ty/D");
	tree->Branch("ppac_xflag", &ppac_xflag, "pxflag/s");
	tree->Branch("ppac_yflag", &ppac_yflag, "pyflag/s");
	tree->Branch("ppac_x", ppac_x, "ppacx[3]/D");
	tree->Branch("ppac_y", ppac_y, "ppacy[3]/D");
	tree->Branch("run", &run, "run/I");
}

}	// ribll
