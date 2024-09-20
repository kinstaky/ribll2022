#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <Math/Vector3D.h>

#include "include/event/threebody_info_event.h"
#include "include/event/generate_event.h"
#include "include/calculator/d2_energy_calculator.h"
#include "include/calculator/target_energy_calculator.h"
#include "include/ppac_track.h"

using namespace ribll;

constexpr double be_excited_energy[4] = {
	0.0, 3.368, 6.179, 7.542
};

int main() {
	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody-sim-0002.root", kGenerateDataPath, kInformationDir
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
	// input event
	ThreeBodyInfoEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// add generated event as friend
	// generated file name
	TString generate_file_name = TString::Format(
		"%s%sgenerate-0002.root", kGenerateDataPath, kSimulateDir
	);
	// add friend
	ipt->AddFriend("gen=tree", generate_file_name);
	// generate event
	GenerateEvent generate_event;
	// setup input branches
	generate_event.SetupInput(ipt, "gen.");

	// output file name
	TString output_file_name = TString::Format(
		"%s%scalculate-t0d2.root", kGenerateDataPath, kShowDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	TH1F ex_diff_c[4];
	for (int i = 0; i < 4; ++i) {
		ex_diff_c[i] = TH1F(
			TString::Format("hexdc%d", i), "difference of excited energy",
			100, -1, 1
		);
	}
	TH2F pid_d1d2("hpd1d2", "D1-D2 PID", 1000, 0, 400, 1000, 0, 400);
	TH2F pid_d2d3("hpd2d3", "D2-D3 PID", 1000, 0, 400, 1000, 0, 400);
	TH2F pid_d1d2_calc("hpd1d2c", "D1-D2 calculated PID", 1000, 0, 400, 1000, 0, 400);
	TH2F pid_d2d3_calc("hpd2d3c", "D2-D3 calculated PID", 1000, 0, 400, 1000, 0, 400);

	// output tree
	TTree opt("tree", "calculate T0D2 energy");
	// output data
	int ppac_flag, taf_flag, bind, hole, target_flag;
	// kinetic energy of particles, indexes for
	// 0: without calculation
	// 1: only calculate the hole part
	// 2: calculate all events
	// 3: calculated Be and measured He
	double be_kinetic[4], he_kinetic[4], d_kinetic, c_kinetic[4];
	// momentum value
	double c_momentum[4];
	// Q value
	double q[4];
	// kinetic energy with target energy lost
	double be_kinetic_target[4], he_kinetic_target[4], d_kinetic_target;
	double c_kinetic_target[4];
	// excited energy ignore 10Be state
	double stateless_excited_energy[3][4];
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
	opt.Branch("c_momentum", c_momentum, "cp[4]/D");
	opt.Branch("q", q, "q[4]/D");
	opt.Branch("be_kinetic_target", be_kinetic_target, "bekt[4]/D");
	opt.Branch("he_kinetic_target", he_kinetic_target, "hekt[4]/D");
	opt.Branch("d_kinetic_target", &d_kinetic_target, "dkt/D");
	opt.Branch("c_kinetic_target", c_kinetic_target, "ckt[4]/D");
	opt.Branch(
		"stateless_excited_energy", stateless_excited_energy, "slext[3][4]/D"
	);


	// T0D2 energy calculators
	elc::D2EnergyCalculator be_calculator("10Be");
	elc::D2EnergyCalculator he_calculator("4He");
	// target energy calculator
	elc::TargetEnergyCalculator be10_target("10Be", "CD2", 9.53);
	elc::TargetEnergyCalculator he4_target("4He", "CD2", 9.53);
	elc::TargetEnergyCalculator h2_target("2H", "CD2", 9.53);

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Showing   0%%");
	fflush(stdout);
	// loop to process
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get data
		ipt->GetEntry(entry);

		// get flags
		taf_flag = event.taf_flag;
		ppac_flag = event.ppac_flag;
		bind = event.bind;
		hole = int(event.hole[0]) + (int(event.hole[1]) << 1);
		target_flag = event.target_flag;

		// ignore special events
		if (taf_flag != 0 || (ppac_flag & 1) == 0) {
			opt.Fill();
			continue;
		}

		// target position
		double tx = event.xptx;
		double ty = event.xpty;

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
		// sum up
		be_kinetic[0] = be_d1_energy + be_d2_energy + be_d3_energy;
		// consider energy loss in target
		be_kinetic_target[0] =
			be10_target.Energy(-0.5/cos(d_be.Theta()), be_kinetic[0]);


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
		// 4He SSD energy
		double he_ssd_energy = 0.0;
		for (int i = 0; i < 3; ++i) {
			// check layer
			if (event.layer[1] <= i+2) break;
			he_ssd_energy +=
				t0_param[i+3][0] + t0_param[i+3][1] * event.ssd_channel[i];
		}
		// sum up
		// without calculation
		he_kinetic[0] =
			he_d1_energy + he_d2_energy + he_d3_energy + he_ssd_energy;

		// deutron kinetic energy
		d_kinetic = event.taf_energy;
		// consider energy lost in target
		d_kinetic_target = h2_target.Energy(-0.5/cos(d_d.Theta()), d_kinetic);



		// calculated 10Be T0D2 energy
		double calc_be_d2_energy = event.layer[0] == 1
			? be_calculator.Energy(0, be_d1_energy, d_be.Theta(), true)
			: be_calculator.DeltaEnergy(1, be_d3_energy, d_be.Theta(), true);
		// must with calculation
		be_kinetic[2] = be_d1_energy + calc_be_d2_energy + be_d3_energy;

		// 4He stop in T0D2
		bool he_stop_d2 = event.layer[1] == 2;
		// 4He calculated T0D2 energy
		double calc_he_d2_energy = event.layer[1] == 1
			? he_calculator.Energy(0, he_d1_energy, d_he.Theta(), true)
			: he_calculator.DeltaEnergy(1, he_d3_energy, d_he.Theta(), he_stop_d2);

		// must with calculation
		he_kinetic[2] =
			he_d1_energy + calc_he_d2_energy + he_d3_energy + he_ssd_energy;
		// consider energy lost in target
		he_kinetic_target[0] =
			he4_target.Energy(-0.5/cos(d_he.Theta()), he_kinetic[0]);


		// calculate energy
		double using_ppac_xz[3] = {ppac_xz[0], ppac_xz[1], ppac_xz[2]};
		double using_ppac_yz[3] = {ppac_yz[0], ppac_yz[1], ppac_yz[2]};
		if (event.run >= ppac_change_run) {
			using_ppac_xz[0] = all_ppac_xz[1];
			using_ppac_yz[0] = all_ppac_yz[1];
		}
		double calc_tx, calc_ty;
		for (int i = 0; i < 2; ++i) {
			// calculate iteration reaction point
			TrackPpac(
				event.xppac_xflag, using_ppac_xz, event.xppac_x,
				be_kinetic[2], he_kinetic[2], d_kinetic,
				event.be_x[0], event.he_x[0], event.d_x, event.d_y,
				calc_tx
			);
			TrackPpac(
				event.xppac_yflag, using_ppac_yz, event.xppac_y,
				be_kinetic[2], he_kinetic[2], d_kinetic,
				event.be_y[0], event.he_y[0], event.d_y, event.d_x,
				calc_ty
			);

			// calculate again
			// 10Be direction vector
			d_be = ROOT::Math::XYZVector(
				event.be_x[0] - calc_tx,
				event.be_y[0] - calc_ty,
				100.0
			).Unit();
			// 4He direction vector
			d_he = ROOT::Math::XYZVector(
				event.he_x[0] - calc_tx,
				event.he_y[0] - calc_ty,
				100.0
			).Unit();
			// 2H direction vector
			d_d = ROOT::Math::XYZVector(
				event.d_x - calc_tx,
				event.d_y - calc_ty,
				135.0
			).Unit();

			// calculated 10Be T0D2 energy
			calc_be_d2_energy = event.layer[0] == 1
				? be_calculator.Energy(0, be_d1_energy, d_be.Theta(), true)
				: be_calculator.DeltaEnergy(1, be_d3_energy, d_be.Theta(), true);
			// must with calculation
			be_kinetic[2] = be_d1_energy + calc_be_d2_energy + be_d3_energy;
			// consider energy lost in target
			be_kinetic_target[2] =
				be10_target.Energy(-0.5/cos(d_be.Theta()), be_kinetic[2]);
			// 4He stop in T0D2?
			he_stop_d2 = event.layer[1] == 2;
			// 4He calculated T0D2 energy
			calc_he_d2_energy = event.layer[1] == 1
				? he_calculator.Energy(0, he_d1_energy, d_he.Theta(), true)
				: he_calculator.DeltaEnergy(
					1, he_d3_energy, d_he.Theta(), he_stop_d2
				);
			// must with calculation
			he_kinetic[2] =
				he_d1_energy + calc_he_d2_energy + he_d3_energy + he_ssd_energy;
			// consider energy lost in target
			he_kinetic_target[2] =
				he4_target.Energy(-0.5/cos(d_he.Theta()), he_kinetic[2]);
		}

		// get target flag
		target_flag = 0;
		if (pow(calc_tx, 2.0) + pow(calc_ty, 2.0) < 200.0) target_flag |= 1;

		// with calculation if hole
		be_kinetic[1] = event.hole[0] ? be_kinetic[2] : be_kinetic[0];
		be_kinetic[3] = be_kinetic[2];
		// He kinetic energy with calculation if needed (hole)
		he_kinetic[1] = event.hole[1] ? he_kinetic[2] : he_kinetic[0];
		he_kinetic[3] = he_kinetic[0];


		// considering energy lost in target
		be_kinetic_target[1] =
			event.hole[0] ? be_kinetic_target[2] : be_kinetic_target[0];
		be_kinetic_target[3] = be_kinetic_target[2];

		he_kinetic_target[1] =
			event.hole[1] ? he_kinetic_target[2] : he_kinetic_target[0];
		he_kinetic_target[3] = he_kinetic_target[0];

		// fill PID
		pid_d1d2.Fill(be_d2_energy, be_d1_energy);
		pid_d1d2_calc.Fill(calc_be_d2_energy, be_d1_energy);
		pid_d1d2.Fill(he_d2_energy, he_d1_energy);
		pid_d1d2_calc.Fill(calc_he_d2_energy, he_d1_energy);
		if (event.layer[0] >= 2) {
			pid_d2d3.Fill(be_d3_energy, be_d2_energy);
			pid_d2d3_calc.Fill(be_d3_energy, calc_be_d2_energy);
		}
		if (event.layer[1] >= 2) {
			pid_d2d3.Fill(he_d3_energy, he_d2_energy);
			pid_d2d3_calc.Fill(he_d3_energy, he_d2_energy);
		}


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
			c_momentum[i] = p_c.R();
			// 14C kinematic energy
			c_kinetic[i] =
				sqrt(pow(c_momentum[i], 2.0) + pow(mass_14c, 2.0)) - mass_14c;

			// three body Q value
			q[i] = be_kinetic[i] + he_kinetic[i] + d_kinetic - c_kinetic[i];

			// ignore state
			for (int j = 0; j < 3; ++j) {
				double stateless_be10_ex = be_excited_energy[j];
				// 10Be momentum
				double be_momentum_stateless = MomentumFromKinetic(
					mass_10be+stateless_be10_ex, be_kinetic_target[i]
				);
				ROOT::Math::XYZVector p_be_stateless =
					d_be * be_momentum_stateless;
				// 4He momentum
				double he_momentum_stateless = MomentumFromKinetic(
					mass_4he, he_kinetic_target[i]
				);
				ROOT::Math::XYZVector p_he_stateless =
					d_he * he_momentum_stateless;

				// excited 14C momentum vector
				ROOT::Math::XYZVector p_excited_c_stateless =
					p_be_stateless + p_he_stateless;
				// excited 14C momentum
				double excited_c_momentum_stateless =
					p_excited_c_stateless.R();

				// excited 14C total energy
				double excited_c_energy_stateless =
					(be_kinetic_target[i] + mass_10be + stateless_be10_ex)
					+ (he_kinetic_target[i] + mass_4he);
				// excited 14C mass
				double excited_c_mass_stateless = sqrt(
					pow(excited_c_energy_stateless, 2.0)
					- pow(excited_c_momentum_stateless, 2.0)
				);
				// excited energy of 14C
				stateless_excited_energy[j][i] =
					excited_c_mass_stateless - mass_14c;
			}
		}

		for (int i = 0; i < 4; ++i) {
			if (q[i] < -10 && q[i] > -14) {
				ex_diff_c[i].Fill(
					stateless_excited_energy[0][i]
					- generate_event.beam_excited_energy
				);
			} else if (q[i] < -14.5 && q[i] > -16) {
				ex_diff_c[i].Fill(
					stateless_excited_energy[1][i]
					- generate_event.beam_excited_energy
				);
			} else if (q[i] < -17.5 && q[i] > -20.0) {
				ex_diff_c[i].Fill(
					stateless_excited_energy[2][i]
					- generate_event.beam_excited_energy
				);
			}
		}

		// fill to tree
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save
	opf.cd();
	pid_d1d2.Write();
	pid_d2d3.Write();
	pid_d1d2_calc.Write();
	pid_d2d3_calc.Write();
	ex_diff_c[0].Write();
	ex_diff_c[2].Write();
	ex_diff_c[3].Write();
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}

