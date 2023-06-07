#include "include/event/channel_event.h"

namespace ribll {

void ChannelEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	// fragments
	tree->SetBranchAddress((prefix+"num").c_str(), &num);
	tree->SetBranchAddress((prefix+"charge").c_str(), daughter_charge);
	tree->SetBranchAddress((prefix+"mass").c_str(), daughter_mass);
	tree->SetBranchAddress((prefix+"energy").c_str(), daughter_energy);
	tree->SetBranchAddress((prefix+"time").c_str(), daughter_time);
	tree->SetBranchAddress((prefix+"px").c_str(), daughter_px);
	tree->SetBranchAddress((prefix+"py").c_str(), daughter_py);
	tree->SetBranchAddress((prefix+"pz").c_str(), daughter_pz);
	tree->SetBranchAddress((prefix+"r").c_str(), daughter_r);
	tree->SetBranchAddress((prefix+"theta").c_str(), daughter_theta);
	tree->SetBranchAddress((prefix+"phi").c_str(), daughter_phi);
	// beam
	tree->SetBranchAddress((prefix+"benergy").c_str(), &beam_energy);
	tree->SetBranchAddress((prefix+"btime").c_str(), &beam_time);
	tree->SetBranchAddress((prefix+"bpx").c_str(), &beam_px);
	tree->SetBranchAddress((prefix+"bpy").c_str(), &beam_py);
	tree->SetBranchAddress((prefix+"bpz").c_str(), &beam_pz);
	tree->SetBranchAddress((prefix+"br").c_str(), &beam_r);
	tree->SetBranchAddress((prefix+"btheta").c_str(), &beam_theta);
	tree->SetBranchAddress((prefix+"bphi").c_str(), &beam_phi);
	// parent
	tree->SetBranchAddress((prefix+"pcharge").c_str(), &parent_charge);
	tree->SetBranchAddress((prefix+"pmass").c_str(), &parent_mass);
	tree->SetBranchAddress((prefix+"penergy").c_str(), &parent_energy);
	tree->SetBranchAddress((prefix+"ppx").c_str(), &parent_px);
	tree->SetBranchAddress((prefix+"ppy").c_str(), &parent_py);
	tree->SetBranchAddress((prefix+"ppz").c_str(), &parent_pz);
	tree->SetBranchAddress((prefix+"pr").c_str(), &parent_r);
	tree->SetBranchAddress((prefix+"ptheta").c_str(), &parent_theta);
	tree->SetBranchAddress((prefix+"pphi").c_str(), &parent_phi);
	// recoil
	tree->SetBranchAddress((prefix+"recoil").c_str(), &recoil);
	tree->SetBranchAddress((prefix+"rcharge").c_str(), &recoil_charge);
	tree->SetBranchAddress((prefix+"rmass").c_str(), &recoil_mass);
	tree->SetBranchAddress((prefix+"renergy").c_str(), &recoil_energy);
	tree->SetBranchAddress((prefix+"rtime").c_str(), &recoil_time);
	tree->SetBranchAddress((prefix+"rpx").c_str(), &recoil_px);
	tree->SetBranchAddress((prefix+"rpy").c_str(), &recoil_py);
	tree->SetBranchAddress((prefix+"rpz").c_str(), &recoil_pz);
	tree->SetBranchAddress((prefix+"rr").c_str(), &recoil_r);
	tree->SetBranchAddress((prefix+"rtheta").c_str(), &recoil_theta);
	tree->SetBranchAddress((prefix+"rphi").c_str(), &recoil_phi);
	// entry
	tree->SetBranchAddress((prefix+"entry").c_str(), &entry);
	tree->SetBranchAddress((prefix+"taf_index").c_str(), &taf_index);
}


void ChannelEvent::SetupOutput(TTree *tree) {
	// fragments
	tree->Branch("num", &num, "num/s");
	tree->Branch("charge", daughter_charge, "Z[num]/s");
	tree->Branch("mass", daughter_mass, "A[num]/s");
	tree->Branch("energy", daughter_energy, "e[num]/D");
	tree->Branch("time" , daughter_time, "t[num]/D");
	tree->Branch("px", daughter_px, "px[num]/D");
	tree->Branch("py", daughter_py, "py[num]/D");
	tree->Branch("pz", daughter_pz, "pz[num]/D");
	tree->Branch("r", daughter_r, "r[num]/D");
	tree->Branch("theta", daughter_theta, "theta[num]/D");
	tree->Branch("phi", daughter_phi, "phi[num]/D");
	// beam
	tree->Branch("benergy", &beam_energy, "be/D");
	tree->Branch("btime", &beam_time, "bt/D");
	tree->Branch("bpx", &beam_px, "bpx/D");
	tree->Branch("bpy", &beam_py, "bpy/D");
	tree->Branch("bpz", &beam_pz, "bpz/D");
	tree->Branch("br", &beam_r, "br/D");
	tree->Branch("btheta", &beam_theta, "btheta/D");
	tree->Branch("bphi", &beam_phi, "bphi/D");
	// parent
	tree->Branch("pcharge", &parent_charge, "pZ/s");
	tree->Branch("pmass", &parent_mass, "pA/s");
	tree->Branch("penergy", &parent_energy, "pe/D");
	tree->Branch("ppx", &parent_px, "ppx/D");
	tree->Branch("ppy", &parent_py, "ppy/D");
	tree->Branch("ppz", &parent_pz, "ppz/D");
	tree->Branch("pr", &parent_r, "pr/D");
	tree->Branch("ptheta", &parent_theta, "ptheta/D");
	tree->Branch("pphi", &parent_phi, "pphi/D");
	// recoil
	tree->Branch("recoil", &recoil, "recoil/s");
	tree->Branch("rcharge", &recoil_charge, "rZ/s");
	tree->Branch("rmass", &recoil_mass, "rA/s");
	tree->Branch("renergy", &recoil_energy, "re/D");
	tree->Branch("rtime", &recoil_time, "rt/D");
	tree->Branch("rpx", &recoil_px, "rpx/D");
	tree->Branch("rpy", &recoil_py, "rpy/D");
	tree->Branch("rpz", &recoil_pz, "rpz/D");
	tree->Branch("rr", &recoil_r, "rr/D");
	tree->Branch("rtheta", &recoil_theta, "rtheta/D");
	tree->Branch("rphi", &recoil_phi, "rphi/D");
	// entry
	tree->Branch("entry", &entry, "entry/L");
	tree->Branch("taf_index", &taf_index, "taf_index/I");
}

}	// namespace ribll