#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/threebody_info_event.h"

using namespace ribll;


inline double VelocityFromMomentum(double momentum, double mass) {
	return momentum / sqrt(pow(momentum, 2.0) + pow(mass, 2.0));
}

inline double GammaFromMomentum(double momentum, double mass) {
	return sqrt(1.0 + pow(momentum/mass, 2.0));
}


void CenterMomentum(
	ROOT::Math::XYZVector momentum1,
	ROOT::Math::XYZVector momentum2,
	double mass1,
	double mass2,
	ROOT::Math::XYZVector &momentum1_center,
	ROOT::Math::XYZVector &momentum2_center
) {
	// gamma
	double gamma1 = GammaFromMomentum(momentum1.R(), mass1);
	double gamma2 = GammaFromMomentum(momentum2.R(), mass2);
	// total energy
	double total_energy1 = gamma1 * mass1;
	double total_energy2 = gamma2 * mass2;


	// mass of center of mass
	double effect_center_mass = gamma1 * mass1 + gamma2 * mass2;
	// center of mass velocity
	ROOT::Math::XYZVector center_velocity =
		(momentum1 + momentum2) / effect_center_mass;	
	// gamma of center of mass
	double center_gamma = 1.0 / sqrt(1.0 - pow(center_velocity.R(), 2.0));

	// momentum1 project to center of mass momentum
	ROOT::Math::XYZVector momentum1_parallel =
		center_velocity.Unit().Dot(momentum1) * center_velocity.Unit();
	// momentum1 orthometric to center of mass momentum
	ROOT::Math::XYZVector momentum1_ortho = momentum1 - momentum1_parallel;

	// momentum2 project to center of mass momentum
	ROOT::Math::XYZVector momentum2_parallel =
		center_velocity.Unit().Dot(momentum2) * center_velocity.Unit();
	// momentum2 orthometric to center of mass momentum
	ROOT::Math::XYZVector momentum2_ortho = momentum2 - momentum2_parallel;
	
	// transmission to center of mass frame
	ROOT::Math::XYZVector momentum1_center_parallel =
		(
			center_gamma * momentum1_parallel.R()
			- sqrt(pow(center_gamma, 2.0) - 1.0) * total_energy1
		) * momentum1_parallel.Unit();
	ROOT::Math::XYZVector momentum2_center_parallel =
		(
			center_gamma * momentum2_parallel.R()
			- sqrt(pow(center_gamma, 2.0) - 1.0) * total_energy2
		) * momentum2_parallel.Unit();

	// reconstruct momentum in center of mass
	momentum1_center = momentum1_center_parallel + momentum1_ortho;
	momentum2_center = momentum2_center_parallel + momentum2_ortho;
}


int main(int argc, char **argv) {
	if (argc > 2) {
		std::cout << "Usage: " << argv[0] << " [tag]\n"
			<< "  tag    9Be or 10Be\n";
		return -1;
	}
	std::string tag = "10Be";
	if (argc == 2) tag = std::string(argv[1]);

	// input file name
	TString info_file_name = TString::Format(
		"%s%sthreebody-%s.root",
		kGenerateDataPath, kInformationDir, tag.c_str()
	);
	// input file
	TFile ipf(info_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< info_file_name << " failed.\n";
		return -1;
	}
	// input event
	ThreeBodyInfoEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// spectrum V2 file name
	TString spectrum_v2_file_name = TString::Format(
		"%s%sthreebody-%s-2.root", kGenerateDataPath, kSpectrumDir, tag.c_str()
	);
	// add friend
	ipt->AddFriend("s=tree", spectrum_v2_file_name);
	// input data
	double spectrum_c_kinetic[4], spectrum_c_momentum[4];
	double spectrum_q[4];
	int spectrum_be_state[4];
	double spectrum_excited_energy[4];
	// setup input branches
	ipt->SetBranchAddress("s.c_kinetic", spectrum_c_kinetic);
	ipt->SetBranchAddress("s.c_momentum", spectrum_c_momentum);
	ipt->SetBranchAddress("s.q", spectrum_q);
	ipt->SetBranchAddress("s.be_state", spectrum_be_state);
	ipt->SetBranchAddress("s.excited_energy_target", spectrum_excited_energy);

	// output file name
	TString output_file_name = TString::Format(
		"%s%sthreebody-extra-%s.root",
		kGenerateDataPath, kInformationDir, tag.c_str()
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "extra threebody information");
	// output data
	bool valid;
	int hole;
	int layer[2];
	double cd_angle;
	double cd_velocity;
	double cd_center_velocity;
	bool in_target;
	bool in_center;
	double c14_kinetic_sigma;
	double c14_momentum_sigma;
	double q, c14_excited_energy;
	int be_state;
	// setup output branches
	opt.Branch("valid", &valid, "valid/O");
	opt.Branch("hole", &hole, "hole/I");
	opt.Branch("layer", layer, "layer[2]/I");
	opt.Branch("cd_angle", &cd_angle, "cda/D");
	opt.Branch("cd_velocity", &cd_velocity, "cdv/D");
	opt.Branch("cd_center_velocity", &cd_center_velocity, "cdcv/D");	
	opt.Branch("in_target", &in_target, "it/O");
	opt.Branch("in_center", &in_center, "ic/O");
	opt.Branch("c14_kinetic_sigma", &c14_kinetic_sigma, "c14ks/D");
	opt.Branch("c14_momentum_sigma", &c14_momentum_sigma, "c14ps/D");
	opt.Branch("q", &q, "q/D");
	opt.Branch("be_state", &be_state, "bes/I");
	opt.Branch("c14_excited_energy", &c14_excited_energy, "c14ex/D");

	// loop
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);
		// initialize
		valid = true;
		hole = 0;
		hole |= event.hole[0] ? 1 : 0;
		hole |= event.hole[1] ? 2 : 0;
		layer[0] = event.layer[0];
		layer[1] = event.layer[1];
		in_target = false;
		in_center = false;
		c14_kinetic_sigma = 10.0;
		c14_momentum_sigma = 10.0;

		// check PPAC
		if (
			(event.ppac_flag & 1) == 0
			|| event.xppac_track[0] == 0
			|| event.xppac_track[1] == 0
		) valid = false;
		// check TAF flag
		if (event.taf_flag != 0) valid = false;
		// check binding events
		if (event.bind != 0) valid = false;
		
		if (!valid) {
			opt.Fill();
			continue;
		}

		// 10Be momentum
		double be_momentum = MomentumFromKinetic(mass_10be, event.t0_energy[0]);
		// 10Be momentum vector direction
		ROOT::Math::XYZVector d_be(
			event.be_x[0] - event.xptx,
			event.be_y[0] - event.xpty,
			100.0
		);
		d_be = d_be.Unit();
		// 10Be momentum vector
		ROOT::Math::XYZVector p_be = d_be * be_momentum;

		// 4He momentum
		double he_momentum = MomentumFromKinetic(mass_4he, event.t0_energy[1]);
		// 4He momentum vector direction
		ROOT::Math::XYZVector d_he(
			event.he_x[0] - event.xptx,
			event.he_y[0] - event.xpty,
			100.0
		);
		d_he = d_he.Unit();
		// 4He momentum vector
		ROOT::Math::XYZVector p_he = d_he * he_momentum;

		// 2H momentum
		double d_momentum = MomentumFromKinetic(mass_2h, event.taf_energy);
		// 2H momentum vector
		ROOT::Math::XYZVector d_d(
			event.d_x - event.xptx,
			event.d_y - event.xpty,
			135.0
		);
		d_d = d_d.Unit();
		// 2H momoentum vector
		ROOT::Math::XYZVector p_d = d_d * d_momentum;

		// excited 14C momentum vector
		ROOT::Math::XYZVector p_excited_c = p_be + p_he;
		// excited 14C momentum
		double excited_c_momentum = p_excited_c.R();
		
		// 14C and 2H angle in lab frame
		cd_angle = acos(p_excited_c.Unit().Dot(p_d.Unit()));

		// excited 14C velocity value

		double velocity_c = VelocityFromMomentum(excited_c_momentum, mass_14c);
		// excited 14C velocity vector
		ROOT::Math::XYZVector v_c = p_excited_c.Unit() * velocity_c;
		
		// 2H velocity value
		double velocity_d = VelocityFromMomentum(d_momentum, mass_2h);
		// 2H velocity vector
		ROOT::Math::XYZVector v_d = p_d.Unit() * velocity_d;

		// relative velocity vector
		ROOT::Math::XYZVector cd_velocity_vec = v_c - v_d;
		cd_velocity = cd_velocity_vec.R();

		// transfrom to center of mass coordinate
		ROOT::Math::XYZVector pcc, pdc;
		CenterMomentum(p_excited_c, p_d, mass_14c, mass_2h, pcc, pdc);

		// get velocity
		ROOT::Math::XYZVector vcc = pcc.Unit() * VelocityFromMomentum(pcc.R(), mass_14c);
		ROOT::Math::XYZVector vdc = pdc.Unit() * VelocityFromMomentum(pdc.R(), mass_2h);

		// relative velocity
		cd_center_velocity = (vcc - vdc).R();

		// check reaction position
		in_center =
			event.xptx > -10 && event.xptx < -5
			&& event.xpty > -3 && event.xpty < 3;
		
		// get kinetic energy and momentum
		double c14_kinetic = spectrum_c_kinetic[3];
		double c14_momentum = spectrum_c_momentum[3];
		if (in_center) {
			c14_kinetic_sigma = (c14_kinetic - 376.0) / 7.4;
			c14_momentum_sigma = (c14_momentum - 3147.0) / 37.3;
		} else {
			c14_kinetic_sigma = (c14_kinetic - 384.0) / 4.6;
			c14_momentum_sigma = (c14_momentum - 3184.0) / 23.7;
		}
		q = spectrum_q[3];
		be_state = spectrum_be_state[3];
		c14_excited_energy = spectrum_excited_energy[3];

		// fill
		opt.Fill();
	}

	// save
	opf.cd();
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();

	return 0;
}