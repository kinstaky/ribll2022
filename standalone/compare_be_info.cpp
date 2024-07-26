#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/threebody_info_event.h"
#include "include/calculator/d2_energy_calculator.h"

using namespace ribll;

// striaght parameters
constexpr double straight_parameters[] = {
	0.47, -0.042,	// D1D2: A, B
	0.59, -0.042	// D2D3: A, B
};

// Be mean and sigma
constexpr double be_ef_parameters[] = {
	19796.1, 250.0,		// D1D2-9Be: mean, sigma
	20757.5, 214.2,		// D1D2-10Be: mean, sigma
	23569.6, 452.53,	// D2D3-9Be: mean, sigma
	24557.7, 293.941	// D2D3-10Be: mean, sigma
};



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


int main() {
	// energy calculator
	elc::D2EnergyCalculator be10_calculator("10Be");

	// output file name
	TString output_file_name = TString::Format(
		"%s%scompare-Be.root",
		kGenerateDataPath, kInformationDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "compare 9Be and 10Be information");
	// output data
	int mass;
	bool hole;
	bool in_target, in_center;
	double c14_kinetic, c14_momentum;
	double c14_kinetic_sigma, c14_momentum_sigma;
	double be_ef, be_calc_ef;
	double be9_sigma, be9_calc_sigma;
	double be10_sigma, be10_calc_sigma;
	double cd_angle, cd_velocity, cd_center_velocity;
	double q, c14_excited_energy;
	int be_state;
	// setup output branches
	opt.Branch("mass", &mass, "A/I");
	opt.Branch("hole", &hole, "hole/O");
	opt.Branch("in_target", &in_target, "it/O");
	opt.Branch("in_center", &in_center, "ic/O");
	opt.Branch("c14_kinetic", &c14_kinetic, "c14k/D");
	opt.Branch("c14_kinetic_sigma", &c14_kinetic_sigma, "c14ks/D");
	opt.Branch("c14_momentum", &c14_momentum, "c14p/D");
	opt.Branch("c14_momentum_sigma", &c14_momentum_sigma, "c14ps/D");
	opt.Branch("be_ef", &be_ef, "beef/D");
	opt.Branch("be_calc_ef", &be_calc_ef, "becef/D");
	opt.Branch("be9_sigma", &be9_sigma, "be9s/D");
	opt.Branch("be9_calc_sigma", &be9_calc_sigma, "be9cs/D");
	opt.Branch("be10_sigma", &be10_sigma, "be10s/D");
	opt.Branch("be10_calc_sigma", &be10_calc_sigma, "be10cs/D");
	opt.Branch("cd_angle", &cd_angle, "cda/D");
	opt.Branch("cd_velocity", &cd_velocity, "cdv/D");
	opt.Branch("cd_center_velocity", &cd_center_velocity, "cdcv/D");
	opt.Branch("q", &q, "q/D");
	opt.Branch("be_state", &be_state, "bes/I");
	opt.Branch("c14_excited_energy", &c14_excited_energy, "c14ex/D");


	for (mass = 9; mass <= 10; ++mass) {
		// input Be file name
		TString be_file_name = TString::Format(
			"%s%sthreebody-calc-%dBe.root",
			kGenerateDataPath, kInformationDir, mass
		);
		// input file
		TFile ipf(be_file_name, "read");
		// input tree
		TTree *ipt = (TTree*)ipf.Get("tree");
		if (!ipt) {
			std::cerr << "Error: Get tree from "
				<< be_file_name << " failed.\n";
			return -1;
		}
		// input event
		ThreeBodyInfoEvent event;
		// setup input branches
		event.SetupInput(ipt);

		// add spectrum-v2 friend
		// file name
		TString spectrum_v2_file_name = TString::Format(
			"%s%sthreebody-%dBe-2.root", kGenerateDataPath, kSpectrumDir, mass
		);
		// add friend
		ipt->AddFriend("s=tree", spectrum_v2_file_name);
		// friend tree data
		double spectrum_q[4], spectrum_excited_energy[4];
		int spectrum_be_state[4];
		// setup input branches
		ipt->SetBranchAddress("s.q", spectrum_q);
		ipt->SetBranchAddress("s.be_state", spectrum_be_state);
		ipt->SetBranchAddress("s.excited_energy", spectrum_excited_energy);
		
		for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
			// get data
			ipt->GetEntry(entry);

			// initialize
			hole = false;
			in_target = false;
			in_center = false;
			c14_kinetic = 0.0;
			c14_momentum = 0.0;
			be9_sigma = 10.0;
			be10_sigma = 10.0;

			// ignore complex conditions
			if (event.taf_flag != 0) continue;
			if (event.xppac_track[0] == 0 || event.xppac_track[1] == 0) continue;
			if (event.hole[1]) continue;
			if (event.bind != 0) continue;
			if (event.layer[0] != 1) continue;
			
			// set hole
			hole = event.hole[0];

			// check reaction position
			in_center =
				event.xptx > -10 && event.xptx < -5
				&& event.xpty > -3 && event.xpty < 3;

			// get kinetic energy and momentum
			c14_kinetic = event.c14_kinetic;
			c14_momentum = event.c14_momentum;
			if (in_center) {
				c14_kinetic_sigma = (c14_kinetic - 376.0) / 7.4;
				c14_momentum_sigma = (c14_momentum - 3147.0) / 37.3;
			} else {
				c14_kinetic_sigma = (c14_kinetic - 384.0) / 4.6;
				c14_momentum_sigma = (c14_momentum - 3184.0) / 23.7;
			}


			// get sigma in straight PID
			// 10Be direction
			ROOT::Math::XYZVector d_be(
				event.be_x[0] - event.xptx,
				event.be_y[0] - event.xpty,
				100.0
			);
			// T0D1 measured energy channel
			double d1e = event.be_channel[0];
			// T0D2 measured energy channel
			double d2e = event.be_channel[1];
			// calculated energy in T0D2
			double cd2e = be10_calculator.Energy(
				0, t0_param[0][0] + t0_param[0][1]*d1e, d_be.Theta(), true
			);
			// T0D2 calculated energy channel
			cd2e = (cd2e - t0_param[1][0]) / t0_param[1][1];
			// measured Be Ef
			be_ef =
				sqrt(d1e*d2e + straight_parameters[0]*d1e*d1e)
				+ straight_parameters[1]*d2e;
			// calculated Be Ef
			be_calc_ef =
				sqrt(d1e*cd2e + straight_parameters[0]*d1e*d1e)
				+ straight_parameters[1]*cd2e;

			// measured 9Be Ef sigma
			be9_sigma = (be_ef - be_ef_parameters[0]) / be_ef_parameters[1];
			// calculated 9Be Ef sigma
			be9_calc_sigma =
				(be_calc_ef - be_ef_parameters[0]) / be_ef_parameters[1];
			// measured 10Be Ef sigma
			be10_sigma = (be_ef - be_ef_parameters[2]) / be_ef_parameters[3];
			// calculated 10Be Ef sigma
			be10_calc_sigma =
				(be_calc_ef - be_ef_parameters[2]) / be_ef_parameters[3];



			// 10Be momentum
			double be_momentum = MomentumFromKinetic(mass_10be, event.t0_energy[0]);
			// 10Be direction
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

			q = spectrum_q[3];
			be_state = spectrum_be_state[3];
			c14_excited_energy = spectrum_excited_energy[3];

			// fill
			opt.Fill();
		}

		// close input file
		ipf.Close();
	}

	// save tree
	opf.cd();
	opt.Write();
	// close files
	opf.Close();
	
	return 0;
}