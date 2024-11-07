#include "include/event/channel_v2_event.h"

namespace ribll {

void ChannelV2Event::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"valid").c_str(), &valid);
	tree->SetBranchAddress((prefix+"ppac_valid").c_str(), &ppac_valid);
	tree->SetBranchAddress((prefix+"hole").c_str(), &hole);
	tree->SetBranchAddress((prefix+"t0_valid").c_str(), &t0_valid);
	tree->SetBranchAddress((prefix+"tafd_edge").c_str(), &tafd_edge);
	tree->SetBranchAddress((prefix+"tafcsi_valid").c_str(), &tafcsi_valid);
	tree->SetBranchAddress((prefix+"ppac_xflag").c_str(), &ppac_xflag);
	tree->SetBranchAddress((prefix+"ppac_yflag").c_str(), &ppac_yflag);
	tree->SetBranchAddress((prefix+"ppac_xnum").c_str(), &ppac_xnum);
	tree->SetBranchAddress((prefix+"ppac_ynum").c_str(), &ppac_ynum);
	tree->SetBranchAddress((prefix+"ppac_x").c_str(), ppac_x);
	tree->SetBranchAddress((prefix+"ppac_y").c_str(), ppac_y);
	tree->SetBranchAddress((prefix+"tx").c_str(), &tx);
	tree->SetBranchAddress((prefix+"ty").c_str(), &ty);
	tree->SetBranchAddress((prefix+"beam_charge").c_str(), &beam_charge);
	tree->SetBranchAddress((prefix+"beam_mass").c_str(), &beam_mass);
	tree->SetBranchAddress((prefix+"beam_energy").c_str(), &beam_energy);
	tree->SetBranchAddress((prefix+"beam_kinetic").c_str(), &beam_kinetic);
	tree->SetBranchAddress((prefix+"beam_momentum").c_str(), &beam_momentum);
	tree->SetBranchAddress((prefix+"parent_charge").c_str(), &parent_charge);
	tree->SetBranchAddress((prefix+"parent_mass").c_str(), &parent_mass);
	tree->SetBranchAddress((prefix+"parent_momentum").c_str(), &parent_momentum);
	tree->SetBranchAddress((prefix+"taf_index").c_str(), &taf_index);
	tree->SetBranchAddress((prefix+"csi_index").c_str(), &csi_index);
	tree->SetBranchAddress((prefix+"recoil_charge").c_str(), &recoil_charge);
	tree->SetBranchAddress((prefix+"recoil_mass").c_str(), &recoil_mass);
	tree->SetBranchAddress((prefix+"recoil_x").c_str(), &recoil_x);
	tree->SetBranchAddress((prefix+"recoil_y").c_str(), &recoil_y);
	tree->SetBranchAddress((prefix+"recoil_z").c_str(), &recoil_z);
	tree->SetBranchAddress((prefix+"recoil_energy").c_str(), &recoil_energy);
	tree->SetBranchAddress((prefix+"recoil_kinetic").c_str(), &recoil_kinetic);
	tree->SetBranchAddress((prefix+"recoil_momentum").c_str(), &recoil_momentum);
	tree->SetBranchAddress((prefix+"fragment_num").c_str(), &fragment_num);
	tree->SetBranchAddress((prefix+"fragment_charge").c_str(), fragment_charge);
	tree->SetBranchAddress((prefix+"fragment_mass").c_str(), fragment_mass);
	tree->SetBranchAddress((prefix+"fragment_x").c_str(), fragment_x);
	tree->SetBranchAddress((prefix+"fragment_y").c_str(), fragment_y);
	tree->SetBranchAddress((prefix+"fragment_z").c_str(), fragment_z);
	tree->SetBranchAddress((prefix+"fragment_energy").c_str(), fragment_energy);
	tree->SetBranchAddress((prefix+"fragment_kinetic").c_str(), fragment_kinetic);
	tree->SetBranchAddress((prefix+"fragment_momentum").c_str(), fragment_momentum);
	tree->SetBranchAddress((prefix+"run").c_str(), &run);
	tree->SetBranchAddress((prefix+"entry").c_str(), &entry);
	tree->SetBranchAddress((prefix+"tafd_front_strip").c_str(), &tafd_front_strip);
	tree->SetBranchAddress((prefix+"tafd_energy").c_str(), &tafd_energy);
}


void ChannelV2Event::SetupOutput(TTree *tree) {
	tree->Branch("valid", &valid, "v/I");
	tree->Branch("ppac_valid", &ppac_valid, "pv/O");
	tree->Branch("hole", &hole, "hole/I");
	tree->Branch("t0_valid", &t0_valid, "t0v/O");
	tree->Branch("tafd_edge", &tafd_edge, "tafdedge/O");
	tree->Branch("tafcsi_valid", &tafcsi_valid, "cv/O");
	tree->Branch("ppac_xflag", &ppac_xflag, "pxflag/s");
	tree->Branch("ppac_yflag", &ppac_yflag, "pyflag/s");
	tree->Branch("ppac_xnum", &ppac_xnum, "pxnum/s");
	tree->Branch("ppac_ynum", &ppac_ynum, "pynum/s");
	tree->Branch("ppac_x", ppac_x, "ppacx[3]/D");
	tree->Branch("ppac_y", ppac_y, "ppacy[3]/D");
	tree->Branch("tx", &tx, "tx/D");
	tree->Branch("ty", &ty, "ty/D");
	tree->Branch("beam_charge", &beam_charge, "bZ/s");
	tree->Branch("beam_mass", &beam_mass, "bA/s");
	tree->Branch("beam_energy", &beam_energy, "be/D");
	tree->Branch("beam_kinetic", &beam_kinetic, "bk/D");
	tree->Branch("beam_momentum", &beam_momentum, "bp/D");
	tree->Branch("parent_charge", &parent_charge, "pZ/s");
	tree->Branch("parent_mass", &parent_mass, "pA/s");
	tree->Branch("parent_momentum", &parent_momentum, "pp/D");
	tree->Branch("taf_index", &taf_index, "tafi/I");
	tree->Branch("csi_index", &csi_index, "ci/I");
	tree->Branch("recoil_charge", &recoil_charge, "rZ/s");
	tree->Branch("recoil_mass", &recoil_mass, "rA/s");
	tree->Branch("recoil_x", &recoil_x, "rx/D");
	tree->Branch("recoil_y", &recoil_y, "ry/D");
	tree->Branch("recoil_z", &recoil_z, "rz/D");
	tree->Branch("recoil_energy", &recoil_energy, "re/D");
	tree->Branch("recoil_kinetic", &recoil_kinetic, "rk/D");
	tree->Branch("recoil_momentum", &recoil_momentum, "rp/D");
	tree->Branch("fragment_num", &fragment_num, "fnum/I");
	tree->Branch("fragment_charge", fragment_charge, "fZ[fnum]/s");
	tree->Branch("fragment_mass", fragment_mass, "fA[fnum]/s");
	tree->Branch("fragment_x", fragment_x, "fx[fnum]/D");
	tree->Branch("fragment_y", fragment_y, "fy[fnum]/D");
	tree->Branch("fragment_z", fragment_z, "fy[fnum]/D");
	tree->Branch("fragment_energy", fragment_energy, "fe[fnum]/D");
	tree->Branch("fragment_kinetic", fragment_kinetic, "fk[fnum]/D");
	tree->Branch("fragment_momentum", fragment_momentum, "fp[fnum]/D");
	tree->Branch("run", &run, "run/S");
	tree->Branch("entry", &entry, "entry/L");
	tree->Branch("tafd_front_strip", &tafd_front_strip, "tafdfs/I");
	tree->Branch("tafd_energy", &tafd_energy, "tafde/D");
}

}	// ribll