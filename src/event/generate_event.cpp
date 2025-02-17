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
		(prefix + "excited_fragment0_kinetic_in_target").c_str(),
		&excited_fragment0_kinetic_in_target
	);
	tree->SetBranchAddress(
		(prefix + "excited_fragment0_theta").c_str(),
		&excited_fragment0_theta
	);
	tree->SetBranchAddress(
		(prefix + "excited_fragment0_phi").c_str(),
		&excited_fragment0_phi
	);
	tree->SetBranchAddress(
		(prefix + "excited_fragment0_x").c_str(), &excited_fragment0_x
	);
	tree->SetBranchAddress(
		(prefix + "excited_fragment0_y").c_str(), &excited_fragment0_y
	);
	tree->SetBranchAddress(
		(prefix + "excited_fragment0_z").c_str(), &excited_fragment0_z
	);
	tree->SetBranchAddress(
		(prefix + "photon_energy").c_str(), &photon_energy
	);
	tree->SetBranchAddress(
		(prefix + "photon_theta").c_str(), &photon_theta
	);
	tree->SetBranchAddress(
		(prefix + "photon_phi").c_str(), &photon_phi
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
		(prefix+"breakup_angle").c_str(), &breakup_angle
	);
	tree->SetBranchAddress(
		(prefix+"gamma_angle").c_str(), &gamma_angle
	);
	tree->SetBranchAddress(
		(prefix+"parent_recoil_angle").c_str(), &parent_recoil_angle
	);
	tree->SetBranchAddress(
		(prefix+"fragment_phi_center").c_str(), &fragment_phi_center
	);
	tree->SetBranchAddress(
		(prefix+"fragment_fragment_angle").c_str(), &fragment_fragment_angle
	);
	tree->SetBranchAddress(
        (prefix+"angle_theta_star").c_str(), &angle_theta_star
    );
	tree->SetBranchAddress(
        (prefix+"angle_psi").c_str(), &angle_psi
    );
	tree->SetBranchAddress(
        (prefix+"angle_chi").c_str(), &angle_chi
    );
	tree->SetBranchAddress(
		(prefix+"recoil_vx").c_str(), &recoil_vx
	);
	tree->SetBranchAddress(
        (prefix+"recoil_vy").c_str(), &recoil_vy
    );
	tree->SetBranchAddress(
        (prefix+"recoil_vz").c_str(), &recoil_vz
    );
	tree->SetBranchAddress(
		(prefix+"parent_vx").c_str(), &parent_vx
	);
	tree->SetBranchAddress(
        (prefix+"parent_vy").c_str(), &parent_vy
    );
	tree->SetBranchAddress(
        (prefix+"parent_vz").c_str(), &parent_vz
    );
	tree->SetBranchAddress(
		(prefix+"fragment_vx").c_str(), fragment_vx
	);
	tree->SetBranchAddress(
        (prefix+"fragment_vy").c_str(), fragment_vy
    );
	tree->SetBranchAddress(
        (prefix+"fragment_vz").c_str(), &fragment_vz
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
		"excited_fragment0_kinetic_in_target",
		&excited_fragment0_kinetic_in_target,
		"xfkit/D"
	);
	tree->Branch(
		"excited_fragment0_theta",
		&excited_fragment0_theta,
		"xftheta/D"
	);
	tree->Branch(
		"excited_fragment0_phi",
		&excited_fragment0_phi,
		"xfphi/D"
	);
	tree->Branch("excited_fragment0_x", &excited_fragment0_x, "xfx/D");
	tree->Branch("excited_fragment0_y", &excited_fragment0_y, "xfy/D");
	tree->Branch("excited_fragment0_z", &excited_fragment0_z, "xfz/D");
	tree->Branch("photon_energy", &photon_energy, "ge/D");
	tree->Branch("photon_theta", &photon_theta, "gtheta/D");
	tree->Branch("photon_phi", &photon_phi, "gphi/D");
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
	tree->Branch("breakup_angle", &breakup_angle, "ba/D");
	tree->Branch("gamma_angle", &gamma_angle, "ga/D");
	tree->Branch("parent_recoil_angle", &parent_recoil_angle, "pra/D");
	tree->Branch("fragment_phi_center", &fragment_phi_center, "fpc/D");
	tree->Branch("fragment_fragment_angle", &fragment_fragment_angle, "ffa/D");
	tree->Branch("angle_theta_star", &angle_theta_star, "athetas/D");
	tree->Branch("angle_psi", &angle_psi, "apsi/D");
	tree->Branch("angle_chi", &angle_chi, "achi/D");
	tree->Branch("recoil_vx", &recoil_vx, "rvx/D");
	tree->Branch("recoil_vy", &recoil_vy, "rvy/D");
	tree->Branch("recoil_vz", &recoil_vz, "rvz/D");
	tree->Branch("parent_vx", &parent_vx, "pvx/D");
	tree->Branch("parent_vy", &parent_vy, "pvy/D");
	tree->Branch("parent_vz", &parent_vz, "pvz/D");
	tree->Branch("fragment_vx", fragment_vx, "fvx[2]/D");
	tree->Branch("fragment_vy", fragment_vy, "fvy[2]/D");
	tree->Branch("fragment_vz", fragment_vz, "fvz[2]/D");
}

}	// namespace