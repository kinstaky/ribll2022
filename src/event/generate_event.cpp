#include "include/event/generate_event.h"

namespace ribll {

void GenerateEvent::SetupInput(TTree *tree, const std::string &prefix) {
	tree->SetBranchAddress(
		(prefix+"depth").c_str(), &depth
	);
	tree->SetBranchAddress(
		(prefix+"beam_kinetic_in_target").c_str(),
		&beam_kinetic_in_target
	);
	tree->SetBranchAddress(
		(prefix+"beam_kinetic_before_target").c_str(),
		&beam_kinetic_before_target
	);
	tree->SetBranchAddress(
		(prefix+"c14_excited").c_str(), &beam_excited_energy
	);
	tree->SetBranchAddress(
		(prefix+"be10_excited").c_str(), &fragment_excited_energy
	);
	tree->SetBranchAddress(
		(prefix+"be10_state").c_str(), &fragment_state
	);
	tree->SetBranchAddress(
		(prefix+"beam_theta").c_str(), &beam_theta
	);
	tree->SetBranchAddress(
		(prefix+"beam_phi").c_str(), &beam_phi
	);
	tree->SetBranchAddress(
		(prefix+"target_x").c_str(), &target_x
	);
	tree->SetBranchAddress(
		(prefix+"target_y").c_str(), &target_y
	);
	tree->SetBranchAddress(
		(prefix+"parent_kinetic").c_str(), &parent_kinetic
	);
	tree->SetBranchAddress(
		(prefix+"parent_theta").c_str(), &parent_theta
	);
	tree->SetBranchAddress(
		(prefix+"parent_phi").c_str(), &parent_phi
	);
	tree->SetBranchAddress(
		(prefix+"recoil_kinetic_after_target").c_str(),
		&recoil_kinetic_after_target
	);
	tree->SetBranchAddress(
		(prefix+"recoil_kinetic_in_target").c_str(),
		&recoil_kinetic_in_target
	);
	tree->SetBranchAddress(
		(prefix+"recoil_theta").c_str(), &recoil_theta
	);
	tree->SetBranchAddress(
		(prefix+"recoil_phi").c_str(), &recoil_phi
	);
	tree->SetBranchAddress(
		(prefix+"recoil_x").c_str(), &rx
	);
	tree->SetBranchAddress(
		(prefix+"recoil_y").c_str(), &ry
	);
	tree->SetBranchAddress(
		(prefix+"recoil_z").c_str(), &rz
	);
	tree->SetBranchAddress(
		(prefix+"recoil_r").c_str(), &rr
	);
	tree->SetBranchAddress(
		(prefix+"fragment_kinetic_after_target").c_str(),
		fragment_kinetic_after_target
	);
	tree->SetBranchAddress(
		(prefix+"fragment_kinetic_in_target").c_str(),
		fragment_kinetic_in_target
	);
	tree->SetBranchAddress(
		(prefix+"fragment_theta").c_str(), fragment_theta
	);
	tree->SetBranchAddress(
		(prefix+"fragment_phi").c_str(), fragment_phi
	);
	tree->SetBranchAddress(
		(prefix+"fragment_x").c_str(), fragment_x
	);
	tree->SetBranchAddress(
		(prefix+"fragment_y").c_str(), fragment_y
	);
	tree->SetBranchAddress(
		(prefix+"fragment_z").c_str(), fragment_z
	);
	tree->SetBranchAddress(
		(prefix+"fragment_r").c_str(), fragment_r
	);
	tree->SetBranchAddress(
		(prefix+"elastic_angle").c_str(), &elastic_angle
	);
	tree->SetBranchAddress(
		(prefix+"parent_recoil_angle").c_str(), &parent_recoil_angle
	);
	tree->SetBranchAddress(
		(prefix+"breakup_angle").c_str(), &breakup_angle
	);
	tree->SetBranchAddress(
		(prefix+"fragment_fragment_angle").c_str(), &fragment_fragment_angle
	);
}


void GenerateEvent::SetupOutput(TTree *tree) {
	tree->Branch("depth", &depth, "depth/D");
	tree->Branch(
		"beam_kinetic_in_target", &beam_kinetic_in_target, "bkit/D"
	);
	tree->Branch(
		"beam_kinetic_before_target", &beam_kinetic_before_target, "bkbt/D"
	);
	tree->Branch("c14_excited", &beam_excited_energy, "c14ex/D");
	tree->Branch("be10_excited", &fragment_excited_energy, "be10ex/D");
	tree->Branch("be10_state", &fragment_state, "be10state/I");
	tree->Branch("beam_theta", &beam_theta, "btheta/D");
	tree->Branch("beam_phi", &beam_phi, "bphi/D");
	tree->Branch("target_x", &target_x, "tx/D");
	tree->Branch("target_y", &target_y, "ty/D");
	tree->Branch("parent_kinetic", &parent_kinetic, "pk/D");
	tree->Branch("parent_theta", &parent_theta, "ptheta/D");
	tree->Branch("parent_phi", &parent_phi, "pphi/D");
	tree->Branch(
		"recoil_kinetic_after_target",
		&recoil_kinetic_after_target,
		"rkat/D"
	);
	tree->Branch(
		"recoil_kinetic_in_target",
		&recoil_kinetic_in_target,
		"rkit/D"
	);
	tree->Branch("recoil_theta", &recoil_theta, "rtheta/D");
	tree->Branch("recoil_phi", &recoil_phi, "rphi/D");
	tree->Branch("recoil_x", &rx, "rx/D");
	tree->Branch("recoil_y", &ry, "ry/D");
	tree->Branch("recoil_z", &rz, "rz/D");
	tree->Branch("recoil_r", &rr, "rr/D");
	tree->Branch(
		"fragment_kinetic_after_target",
		fragment_kinetic_after_target,
		"fkat[2]/D"
	);
	tree->Branch(
		"fragment_kinetic_in_target",
		fragment_kinetic_in_target,
		"fkit[2]/D"
	);
	tree->Branch("fragment_theta", fragment_theta, "ftheta[2]/D");
	tree->Branch("fragment_phi", fragment_phi, "fphi[2]/D");
	tree->Branch("fragment_x", fragment_x, "fx[2]/D");
	tree->Branch("fragment_y", fragment_y, "fy[2]/D");
	tree->Branch("fragment_z", fragment_z, "fz[2]/D");
	tree->Branch("fragment_r", fragment_r, "fr[2]/D");
	tree->Branch("elastic_angle", &elastic_angle, "ea/D");
	tree->Branch("parent_recoil_angle", &parent_recoil_angle, "pra/D");
	tree->Branch("breakup_angle", &breakup_angle, "ba/D");
	tree->Branch("fragment_fragment_angle", &fragment_fragment_angle, "ffa/D");
}

}	// namespace