#include "include/event/channel_event.h"

namespace ribll {

void ChannelEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	// fragments
	tree->SetBranchAddress((prefix+"num").c_str(), &num);
	tree->SetBranchAddress((prefix+"charge").c_str(), charge);
	tree->SetBranchAddress((prefix+"mass").c_str(), mass);
	tree->SetBranchAddress((prefix+"energy").c_str(), energy);
	tree->SetBranchAddress((prefix+"px").c_str(), px);
	tree->SetBranchAddress((prefix+"py").c_str(), py);
	tree->SetBranchAddress((prefix+"pz").c_str(), pz);
	tree->SetBranchAddress((prefix+"r").c_str(), r);
	tree->SetBranchAddress((prefix+"theta").c_str(), theta);
	tree->SetBranchAddress((prefix+"phi").c_str(), phi);
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
	tree->SetBranchAddress((prefix+"rpx").c_str(), &recoil_px);
	tree->SetBranchAddress((prefix+"rpy").c_str(), &recoil_py);
	tree->SetBranchAddress((prefix+"rpz").c_str(), &recoil_pz);
	tree->SetBranchAddress((prefix+"rr").c_str(), &recoil_r);
	tree->SetBranchAddress((prefix+"rtheta").c_str(), &recoil_theta);
	tree->SetBranchAddress((prefix+"rphi").c_str(), &recoil_phi);
}


void ChannelEvent::SetupOutput(TTree *tree) {
	// fragments
	tree->Branch("num", &num, "num/s");
	tree->Branch("charge", charge, "Z[num]/s");
	tree->Branch("mass", mass, "A[num]/s");
	tree->Branch("energy", energy, "e[num]/D");
	tree->Branch("px", px, "px[num]/D");
	tree->Branch("py", py, "py[num]/D");
	tree->Branch("pz", pz, "pz[num]/D");
	tree->Branch("r", r, "r[num]/D");
	tree->Branch("theta", theta, "theta[num]/D");
	tree->Branch("phi", phi, "phi[num]/D");
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
	tree->Branch("rcharge", &recoil_charge, "rZ/s");
	tree->Branch("rmass", &recoil_mass, "rA/s");
	tree->Branch("renergy", &recoil_energy, "re/D");
	tree->Branch("rpx", &recoil_px, "rpx/D");
	tree->Branch("rpy", &recoil_py, "rpy/D");
	tree->Branch("rpz", &recoil_pz, "rpz/D");
	tree->Branch("rr", &recoil_r, "rr/D");
	tree->Branch("rtheta", &recoil_theta, "rtheta/D");
	tree->Branch("rphi", &recoil_phi, "rphi/D");
}

}	// namespace ribll