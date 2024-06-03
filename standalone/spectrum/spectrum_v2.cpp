// spectrum 程序版本2，版本1太臃肿了，重新写一个
#include <iostream>
#include <fstream>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <Math/Vector3D.h>

#include "include/event/threebody_info_event.h"
#include "include/calculator/d2_energy_calculator.h"
#include "include/calculator/target_energy_calculator.h"

using namespace ribll;

int main(int argc, char **argv) {
	std::string suffix = "";
	if (argc > 1) {
		suffix = std::string(argv[1]);
	}

	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody%s.root",
		kGenerateDataPath,
		kInformationDir,
		suffix.empty() ? "" : ("-"+suffix).c_str()
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// input data
	ThreeBodyInfoEvent event;
	// setup input branches
	event.SetupInput(ipt);


	// output file name
	TString output_file_name = TString::Format(
		"%s%sthreebody%s-2.root",
		kGenerateDataPath,
		kSpectrumDir,
		suffix.empty() ? "" : ("-"+suffix).c_str()
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "spectrum");
	// output data
	int ppac_flag, taf_flag, bind, hole, target_flag;
	// kinetic energy of particles, indexes for
	// 0: without calculation
	// 1: only calculate the hole part
	// 2: calculate all events
	// 3: calculated Be and measured He
	double be_kinetic[4], he_kinetic[4], d_kinetic, c_kinetic[4];
	// Q value 
	double q[4];
	// 10Be state, under different T0D2 energy method
	int be_state[4];
	// excited energy under different T0D2 energy method
	double excited_energy[4];
	// kinetic energy with target energy lost
	double be_kinetic_target[3], he_kinetic_target[3], d_kinetic_target;
	double c_target_kinetic[3];
	// excited energy with target energy lost
	double excited_energy_target[3];
	// correct the kinetic based on the simulation result
	double be_kinetic_target_correct[3], he_kinetic_target_correct[3];
	double excited_energy_target_correct[3];
	// setup output branches
	opt.Branch("ppac_flag", &ppac_flag, "pflag/I");
	opt.Branch("taf_flag", &taf_flag, "tflag/I");
	opt.Branch("bind", &bind, "bind/I");
	opt.Branch("hole", &hole, "hole/I");
	opt.Branch("target_flag", &target_flag, "tarflag/I");
	opt.Branch("be_kinetic", be_kinetic, "bek[4]/D");
	opt.Branch("he_kinetic", he_kinetic, "hek[4]/D");
	opt.Branch("d_kinetic", &d_kinetic, "dk/D");
	opt.Branch("c_kinetic", c_kinetic, "ck[4]/D");
	opt.Branch("q", q, "q[4]/D");
	opt.Branch("be_state", be_state, "bes[4]/I");
	opt.Branch("excited_energy", excited_energy, "ex[4]/D");
	opt.Branch("be_kinetic_target", be_kinetic_target, "bekt[3]/D");
	opt.Branch("he_kinetic_target", he_kinetic_target, "hekt[3]/D");
	opt.Branch("d_kinetic_target", &d_kinetic_target, "dkt/D");
	opt.Branch("c_target_kinetic", c_target_kinetic, "ckt[3]/D");
	opt.Branch("excited_energy_target", excited_energy_target, "ext[3]/D");
	opt.Branch(
		"be_kinetic_target_correct",
		be_kinetic_target_correct,
		"bektc[3]/D"
	);
	opt.Branch(
		"he_kinetic_target_correct",
		he_kinetic_target_correct,
		"hektc[3]/D"
	);
	opt.Branch(
		"excited_energy_target_correct",
		excited_energy_target_correct,
		"extc[3]/D"
	);

	// Q value correct
	constexpr double q_correct[12] = {
		-0.46, -0.25, 0.18, 0.49,
		0.49, 0.55, 0.28, 0.30,
		-0.04, -0.10, -0.12, -0.64
	};

	// T0D2 energy calculators
	elc::D2EnergyCalculator be_calculator("10Be");
	elc::D2EnergyCalculator he_calculator("4He");
	// target energy calculator
	elc::TargetEnergyCalculator be10_target("10Be", "CD2", 9.53);
	elc::TargetEnergyCalculator he4_target("4He", "CD2", 9.53);
	elc::TargetEnergyCalculator h2_target("2H", "CD2", 9.53);


	// loop to process
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		// get data
		ipt->GetEntry(entry);

		// get flags
		taf_flag = event.taf_flag;
		ppac_flag = event.ppac_flag;
		bind = event.bind;
		hole = int(event.hole[0]) + (int(event.hole[1]) << 1);
		target_flag = event.target_flag;

		// ignore special events
		if (taf_flag != 0) {
			opt.Fill();
			continue;
		}

		// target position
		double tx = (event.ppac_flag & 1) == 0 ? event.vptx : event.xptx;
		double ty = (event.ppac_flag & 1) == 0 ? event.vpty : event.xpty;
		
		// 10Be direction vector
		ROOT::Math::XYZVector d_be(
			event.be_x[0] - tx,
			event.be_y[0] - ty,
			100.0
		);
		d_be = d_be.Unit();
		// 4He direction vector
		ROOT::Math::XYZVector d_he(
			event.he_x[0] - tx,
			event.he_y[0] - ty,
			100.0
		);
		d_he = d_he.Unit();
		// 2H direction vector
		ROOT::Math::XYZVector d_d(
			event.d_x - tx,
			event.d_y - ty,
			135.0
		);
		d_d = d_d.Unit();


		// calculate particle kinetic energy
		// 10Be kinetic energy
		// 10Be T0D1 energy
		double be_d1_energy =
			t0_param[0][0] + t0_param[0][1] * event.be_channel[0];
		// 10Be T0D2 energy
		double be_d2_energy =
			t0_param[1][0] + t0_param[1][1] * event.be_channel[1];
		// 10Be T0D3 energy
		double be_d3_energy = event.layer[0] >= 2
			? t0_param[2][0] + t0_param[2][1] * event.be_channel[2]
			: 0.0;
		// calculated 10Be T0D2 energy
		double calc_be_d2_energy = event.layer[0] == 1
			? be_calculator.Energy(0, be_d1_energy, d_be.Theta(), true)
			: be_calculator.DeltaEnergy(1, be_d3_energy, d_be.Theta(), true);
		// sum up
		// without calculation
		be_kinetic[0] = be_d1_energy + be_d2_energy + be_d3_energy;
		// must with calculation
		be_kinetic[2] = be_d1_energy + calc_be_d2_energy + be_d3_energy;
		// with calculation if hole
		be_kinetic[1] = event.hole[0] ? be_kinetic[2] : be_kinetic[0];
		be_kinetic[3] = be_kinetic[2];

		// 4He kinetic energy
		// 4He T0D1 energy
		double he_d1_energy =
			t0_param[0][0] + t0_param[0][1] * event.he_channel[0];
		// 4He T0D2 energy
		double he_d2_energy =
			t0_param[1][0] + t0_param[1][1] * event.he_channel[1];
		// 4He T0D3 energy
		double he_d3_energy = event.layer[1] >= 2
			? t0_param[2][0] + t0_param[2][1] * event.he_channel[2]
			: 0.0;
		// 4He stop in T0D2
		bool he_stop_d2 = event.layer[1] == 2;
		// 4He calculated T0D2 energy
		double calc_he_d2_energy = event.layer[1] == 1
			? he_calculator.Energy(0, he_d1_energy, d_he.Theta(), true)
			: he_calculator.DeltaEnergy(1, he_d3_energy, d_he.Theta(), he_stop_d2);
		// sum DSSD energy
		// without calculation
		he_kinetic[0] = he_d1_energy + he_d2_energy + he_d3_energy;
		// must with calculation
		he_kinetic[2] = he_d1_energy + calc_he_d2_energy + he_d3_energy;
		// 4He SSD energy
		for (int i = 0; i < 3; ++i) {
			// check layer
			if (event.layer[1] <= i+2) break;
			double ssd_energy =
				t0_param[i+3][0] + t0_param[i+3][1] * event.ssd_channel[i];
			// sum up
			he_kinetic[0] += ssd_energy;
			he_kinetic[2] += ssd_energy;
		}
		// He kinetic energy with calculation if needed (hole)
		he_kinetic[1] = event.hole[1] ? he_kinetic[2] : he_kinetic[0];
		he_kinetic[3] = he_kinetic[0];

		// deutron kinetic energy
		d_kinetic = event.taf_energy;


		// consider energy loss in target
		be_kinetic_target[0] =
			be10_target.Energy(-0.5/cos(d_be.Theta()), be_kinetic[0]);
		be_kinetic_target[2] =
			be10_target.Energy(-0.5/cos(d_be.Theta()), be_kinetic[2]);
		be_kinetic_target[1] =
			event.hole[0] ? be_kinetic_target[2] : be_kinetic_target[0];
		// consider energy lost in target
		he_kinetic_target[0] =
			he4_target.Energy(-0.5/cos(d_he.Theta()), he_kinetic[0]);
		he_kinetic_target[2] =
			he4_target.Energy(-0.5/cos(d_he.Theta()), he_kinetic[2]);
		he_kinetic_target[1] =
			event.hole[1] ? he_kinetic_target[2] : he_kinetic_target[0];
		// consider energy lost in target
		d_kinetic_target = h2_target.Energy(-0.5/cos(d_d.Theta()), d_kinetic);


		for (int i = 0; i < 4; ++i) {
			// 10Be momentum
			double be_momentum = MomentumFromKinetic(mass_10be, be_kinetic[i]);
			ROOT::Math::XYZVector p_be = d_be * be_momentum;

			// 4He momentum
			double he_momentum = MomentumFromKinetic(mass_4he, he_kinetic[i]);
			ROOT::Math::XYZVector p_he = d_he * he_momentum;

			// 2H momentum
			double d_momentum = MomentumFromKinetic(mass_2h, d_kinetic);
			ROOT::Math::XYZVector p_d = d_d * d_momentum;

			// beam 14C momentum vector
			ROOT::Math::XYZVector p_c = p_be + p_he + p_d;

			// 14C momentum
			double c_momentum = p_c.R();
			// 14C kinematic energy
			c_kinetic[i] =
				sqrt(pow(c_momentum, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

			// three body Q value
			q[i] = be_kinetic[i] + he_kinetic[i] + d_kinetic - c_kinetic[i];
			q[i] -= q_correct[event.csi_index];

			// get 10Be state from Q value
			if (q[i] < -11 && q[i] > -13) be_state[i] = 0;
			else if (q[i] < -14.5 && q[i] > -16.2) be_state[i] = 1;
			else if (q[i] < -17.3 && q[i] > -19) be_state[i] = 2;
			// else if (q[i] < -19 && q[i] > -20.5) be_state[i] = 3;
			else be_state[i] = -1;

			double be10_excited_energy = 0.0;
			if (be_state[i] == 1) be10_excited_energy = 3.368;
			else if (be_state[i] == 2) be10_excited_energy = 6.179;
			else if (be_state[i] == 3) be10_excited_energy = 7.542;

			// excited 14C momentum vector
			ROOT::Math::XYZVector p_excited_c = p_be + p_he;
			// excited 14C momentum
			double excited_c_momentum = p_excited_c.R();
			// excited 14C total energy
			double excited_c_energy =
				(be_kinetic[i] + mass_10be + be10_excited_energy)
				+ (he_kinetic[i] + mass_4he);
			// excited 14C mass
			double excited_c_mass = sqrt(
				pow(excited_c_energy, 2.0) - pow(excited_c_momentum, 2.0)
			);
			// excited energy of 14C
			excited_energy[i] = excited_c_mass - mass_14c;
		}

		// consider energy lost in target
		for (int i = 0; i < 3; ++i) {
			// 10Be momentum
			double be_momentum =
				MomentumFromKinetic(mass_10be, be_kinetic_target[i]);
			ROOT::Math::XYZVector p_be = d_be * be_momentum;

			// 4He momentum
			double he_momentum =
				MomentumFromKinetic(mass_4he, he_kinetic_target[i]);
			ROOT::Math::XYZVector p_he = d_he * he_momentum;

			double be10_excited_energy = 0.0;
			if (be_state[i] == 1) be10_excited_energy = 3.368;
			else if (be_state[i] == 2) be10_excited_energy = 6.179;
			else if (be_state[i] == 3) be10_excited_energy = 7.542;

			// excited 14C momentum vector
			ROOT::Math::XYZVector p_excited_c = p_be + p_he;
			// excited 14C momentum
			double excited_c_momentum = p_excited_c.R();
			// excited 14C total energy
			double excited_c_energy =
				(be_kinetic_target[i] + mass_10be + be10_excited_energy)
				+ (he_kinetic_target[i] + mass_4he);
			// excited 14C mass
			double excited_c_mass = sqrt(
				pow(excited_c_energy, 2.0) - pow(excited_c_momentum, 2.0)
			);
			// excited energy of 14C
			excited_energy_target[i] = excited_c_mass - mass_14c;
		}


		// energy lost in target and correct
		if (event.layer[0] == 1) {
			be_kinetic_target_correct[0] = be_kinetic_target[0] - 1.72;
			be_kinetic_target_correct[2] = be_kinetic_target[2] - 0.4;
		} else if (event.layer[0] == 2) {
			be_kinetic_target_correct[0] = be_kinetic_target[0] - 1.39;
			be_kinetic_target_correct[2] = be_kinetic_target[2] - 0.84;
		}
		be_kinetic_target_correct[1] = event.hole[0]
			? be_kinetic_target_correct[2]
			: be_kinetic_target_correct[0];
		if (event.layer[1] == 1) {
			he_kinetic_target_correct[0] = he_kinetic_target[0] - 0.51;
			he_kinetic_target_correct[2] = he_kinetic_target[2] + 0.56;
		} else if (event.layer[1] == 2) {
			he_kinetic_target_correct[0] = he_kinetic_target[0] - 0.39;
			he_kinetic_target_correct[2] = he_kinetic_target[2] + 0.16;
		} else if (event.layer[1] > 2) {
			he_kinetic_target_correct[0] = he_kinetic_target[0] - 0.31;
			he_kinetic_target_correct[2] = he_kinetic_target[2] - 0.20;
		}
		he_kinetic_target_correct[1] = event.hole[1]
			? he_kinetic_target_correct[2]
			: he_kinetic_target_correct[0];
		for (int i = 0; i < 3; ++i) {
			// 10Be momentum
			double be_momentum =
				MomentumFromKinetic(mass_10be, be_kinetic_target_correct[i]);
			ROOT::Math::XYZVector p_be = d_be * be_momentum;

			// 4He momentum
			double he_momentum =
				MomentumFromKinetic(mass_4he, he_kinetic_target_correct[i]);
			ROOT::Math::XYZVector p_he = d_he * he_momentum;

			double be10_excited_energy = 0.0;
			if (be_state[i] == 1) be10_excited_energy = 3.368;
			else if (be_state[i] == 2) be10_excited_energy = 6.179;
			else if (be_state[i] == 3) be10_excited_energy = 7.542;

			// excited 14C momentum vector
			ROOT::Math::XYZVector p_excited_c = p_be + p_he;
			// excited 14C momentum
			double excited_c_momentum = p_excited_c.R();
			// excited 14C total energy
			double excited_c_energy =
				(be_kinetic_target_correct[i]+mass_10be+be10_excited_energy)
				+ (he_kinetic_target_correct[i] + mass_4he);
			// excited 14C mass
			double excited_c_mass = sqrt(
				pow(excited_c_energy, 2.0) - pow(excited_c_momentum, 2.0)
			);
			// excited energy of 14C
			excited_energy_target_correct[i] = excited_c_mass - mass_14c;
		}

		// fill to tree
		opt.Fill();
	}
	
	// save tree
	opf.cd();
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();

	return 0;
}