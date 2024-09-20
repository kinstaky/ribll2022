// spectrum 程序版本2，版本1太臃肿了，重新写一个
#include <iostream>
#include <fstream>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <Math/Vector3D.h>

#include "include/event/threebody_info_event.h"
#include "include/calculator/d2_energy_calculator.h"
#include "include/calculator/target_energy_calculator.h"
#include "include/ppac_track.h"

using namespace ribll;

constexpr double seperate_q_border[12][6] = {
	{-20.0, -17.5, -16.8, -14.8, -13.7, -11.0},	// good
	{-19.4, -17.8, -16.3, -14.5, -13.2, -10.7},	// bad
	{-19.1, -17.3, -16.0, -14.2, -13.2, -10.8},	// good
	{-19.2, -16.8, -15.8, -14.2, -13.0, -11.0}, // not so good
	{-19.0, -16.5, -15.6, -13.4, -12.8, -10.5},	// good
	{-18.7, -16.7, -15.8, -14.0, -12.8, -10.2},	// not so good
	{-19.0, -16.7, -15.8, -14.4, -13.5, -10.5},	// not so good
	{-19.3, -17.4, -15.7, -14.0, -13.4, -10.6}, // not so good
	{-19.6, -17.5, -15.5, -14.2, -13.0, -11.0},	// bad
	{-19.8, -17.6, -16.4, -14.8, -13.8, -11.0}, // low statistics
	{-20.1, -18.4, -16.6, -15.0, -13.8, -11.4},	// bad
	{-19.5, -17.3, -16.7, -14.8, -13.6, -11.2}	// good
};


constexpr double be_excited_energy[4] = {
	0.0, 3.368, 6.179, 7.542
};


// straight parameters
constexpr double sa12 = 0.47;
constexpr double sb12 = -0.042;
constexpr double sa23 = 0.59;
constexpr double sb23 = -0.042;


inline double Straight(double de, double e, double a, double b) {
	return sqrt(de*e + a*de*de) + b*e;
}


int main(int argc, char **argv) {
	std::string suffix = "";
	if (argc > 1) {
		if (argv[1][0] == '-' && argv[1][1] == 'h') {
			std::cout << "Usage: " << argv[0] << " [suffix]\n"
				<< "  suffix        ROOT file suffix tag\n";
			return -1;
		}
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
	// histogram of stateless excited energy
	TH1F hist_stateless_excited_energy[3][11];
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 11; ++j) {
			hist_stateless_excited_energy[i][j] = TH1F(
				TString::Format("hs%dq%d", i, j),
				TString::Format("state %d Q [%d, %d)", i, -21+j, -20+j),
				150, 12, 27
			);
		}
	}
	// PID
	TH2F pid_d1d2("hpd1d2", "D1-D2 PID", 1000, 0, 250, 1000, 0, 200);
	TH2F pid_d2d3("hpd2d3", "D2-D3 PID", 1000, 0, 150, 1000, 0, 250);
	TH2F pid_d1d2_calc("hpd1d2c", "D1-D2 calculated PID", 1000, 0, 250, 1000, 0, 200);
	TH2F pid_d2d3_calc("hpd2d3c", "D2-D3 calculated PID", 1000, 0, 150, 1000, 0, 250);
	// straight PID
	TH1F pid_d1d2_straight("hpd1d2s", "D1-D2 straight PID", 3000, 0, 30000);
	TH1F pid_d2d3_straight("hpd2d3s", "D2-D3 straight PID", 2000, 0, 40000);
	TH1F pid_d1d2_calc_straight("hpd1d2cs", "D1-D2 calculated straight PID", 3000, 0, 30000);
	TH1F pid_d2d3_calc_straight("hpd2d3cs", "D2-D3 calculated straight PID", 2000, 0, 40000);
	// Q value spectrum
	TH1F hq[4];
	for (int i = 0; i < 4; ++i) {
		hq[i] = TH1F(TString::Format("hq%d", i), "Q value", 70, -23, -8);
	}
	// output tree
	TTree opt("tree", "spectrum");
	// output data
	int ppac_flag, taf_flag, bind, hole, target_flag;
	// PID in correct range? index 0 for origin, 1 for calculated 
	int straight[2];
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
	// 10Be state, under different T0D2 energy method
	int be_state[4];
	int sep_be_state[4];
	// excited energy under different T0D2 energy method
	double excited_energy[4];
	// kinetic energy with target energy lost
	double be_kinetic_target[4], he_kinetic_target[4], d_kinetic_target;
	double c_kinetic_target[4];
	// excited energy with target energy lost
	double excited_energy_target[4];
	// 14C excited energy considering target lost and seperated CsI
	double sep_ex_target[4];
	// excited energy ignore 10Be state
	double stateless_excited_energy[3][4];

	// setup output branches
	opt.Branch("ppac_flag", &ppac_flag, "pflag/I");
	opt.Branch("taf_flag", &taf_flag, "tflag/I");
	opt.Branch("bind", &bind, "bind/I");
	opt.Branch("hole", &hole, "hole/I");
	opt.Branch("target_flag", &target_flag, "tarflag/I");
	opt.Branch("straight", &straight, "straight[2]/I");
	opt.Branch("be_kinetic", be_kinetic, "bek[4]/D");
	opt.Branch("he_kinetic", he_kinetic, "hek[4]/D");
	opt.Branch("d_kinetic", &d_kinetic, "dk/D");
	opt.Branch("c_kinetic", c_kinetic, "ck[4]/D");
	opt.Branch("c_momentum", c_momentum, "cp[4]/D");
	opt.Branch("q", q, "q[4]/D");
	opt.Branch("be_state", be_state, "bes[4]/I");
	opt.Branch("sep_be_state", sep_be_state, "sbes[4]/I");
	opt.Branch("excited_energy", excited_energy, "ex[4]/D");
	opt.Branch("be_kinetic_target", be_kinetic_target, "bekt[4]/D");
	opt.Branch("he_kinetic_target", he_kinetic_target, "hekt[4]/D");
	opt.Branch("d_kinetic_target", &d_kinetic_target, "dkt/D");
	opt.Branch("c_kinetic_target", c_kinetic_target, "ckt[4]/D");
	opt.Branch("excited_energy_target", excited_energy_target, "ext[4]/D");
	opt.Branch("spe_ex_target", sep_ex_target, "sext[4]/D");
	opt.Branch(
		"stateless_excited_energy", stateless_excited_energy, "slext[3][4]/D"
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
		straight[0] = straight[1] = 0;

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



		be_kinetic_target[1] =
			event.hole[0] ? be_kinetic_target[2] : be_kinetic_target[0];
		be_kinetic_target[3] = be_kinetic_target[2];

		he_kinetic_target[1] =
			event.hole[1] ? he_kinetic_target[2] : he_kinetic_target[0];
		he_kinetic_target[3] = he_kinetic_target[0];


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

			// get seperated Be state
			sep_be_state[i] = -1;
			for (int j = 0; j < 3; ++j) {
				if (
					q[i] > seperate_q_border[event.csi_index][2*j]
					&& q[i] < seperate_q_border[event.csi_index][2*j+1]
				) {
					sep_be_state[i] = 2-j;
					break;
				}
			}
			double sep_be_ex = sep_be_state[i] >= 0
				? be_excited_energy[sep_be_state[i]]
				: 0.0;


			// correct Q
			q[i] -= q_correct[event.csi_index];
			// get 10Be state from Q value
			if (q[i] < -11 && q[i] > -14) be_state[i] = 0;
			else if (q[i] < -14.5 && q[i] > -16.2) be_state[i] = 1;
			else if (q[i] < -17.3 && q[i] > -19) be_state[i] = 2;
			// else if (q[i] < -19 && q[i] > -20.5) be_state[i] = 3;
			else be_state[i] = -1;
			// get 10Be excited energy form state
			double be10_excited_energy = be_state[i] >= 0
				? be_excited_energy[be_state[i]]
				: 0.0;

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



			// consider energy lost in target
			// 10Be momentum
			double be_momentum_target = MomentumFromKinetic(
				mass_10be+be10_excited_energy, be_kinetic_target[i]
			);
			ROOT::Math::XYZVector p_be_target = d_be * be_momentum_target;

			// 4He momentum
			double he_momentum_target =
				MomentumFromKinetic(mass_4he, he_kinetic_target[i]);
			ROOT::Math::XYZVector p_he_target = d_he * he_momentum_target;


			// excited 14C momentum vector
			ROOT::Math::XYZVector p_excited_c_target = p_be_target + p_he_target;
			// excited 14C momentum
			double excited_c_momentum_target = p_excited_c_target.R();
			// 14C kinematic energy in target
			c_kinetic_target[i] = sqrt(
				pow(excited_c_momentum_target, 2.0) + pow(mass_14c, 2.0)
			) - mass_14c;

			// seperated excited 14C total energy in target
			double sep_ex_c_energy_target =
				(be_kinetic_target[i] + mass_10be + sep_be_ex)
				+ (he_kinetic_target[i] + mass_4he);
			// seperated excited 14C mass
			double sep_ex_c_mass_target = sqrt(
				pow(sep_ex_c_energy_target, 2.0)
				- pow(excited_c_momentum_target, 2.0)
			);
			sep_ex_target[i] = sep_ex_c_mass_target - mass_14c;


			// excited 14C total energy
			double excited_c_energy_target =
				(be_kinetic_target[i] + mass_10be + be10_excited_energy)
				+ (he_kinetic_target[i] + mass_4he);
			// excited 14C mass
			double excited_c_mass_target = sqrt(
				pow(excited_c_energy_target, 2.0)
				- pow(excited_c_momentum_target, 2.0)
			);
			// excited energy of 14C
			excited_energy_target[i] = excited_c_mass_target - mass_14c;

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

		for (int i = 0; i < 3; ++i) {
			if (pow(event.xptx, 2.0) + pow(event.xpty, 2.0) >= 200.0) continue;
			// q index
			int q_index = int(q[3]+21.0);
			if (q_index >= 0 && q_index <= 10) {
				hist_stateless_excited_energy[i][q_index].Fill(
					stateless_excited_energy[i][3]
				);
			}
		}


		// fill histogram to compare different methods
		if (
			(event.ppac_flag & 1) == 1
			&& event.taf_flag == 0
			&& (event.target_flag & 1) == 1
			&& bind == 0
			&& hole == 0
		) {
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
				pid_d2d3_calc.Fill(he_d3_energy, calc_he_d2_energy);
			}

			// fill straight PID
			if (event.layer[0] == 1) {
				double ef = Straight(
					(be_d1_energy - t0_param[0][0]) / t0_param[0][1],
					(be_d2_energy - t0_param[1][0]) / t0_param[1][1],
					sa12, sb12
				);
				pid_d1d2_straight.Fill(ef);
				if (ef > 20200 && ef < 21200) straight[0] |= 1;
				double ef_calc = Straight(
					(be_d1_energy - t0_param[0][0]) / t0_param[0][1],
					(calc_be_d2_energy - t0_param[1][0]) / t0_param[1][1],
					sa12, sb12
				);
				pid_d1d2_calc_straight.Fill(ef_calc);
				if (ef_calc > 20200 && ef_calc < 21200) straight[1] |= 1;
			}
			if (event.layer[1] == 1) {
				double ef = Straight(
					(he_d1_energy - t0_param[0][0]) / t0_param[0][1],
					(he_d2_energy - t0_param[1][0]) / t0_param[1][1],
					sa12, sb12
				);
				pid_d1d2_straight.Fill(ef);
				if (ef > 6100 && ef < 6600) straight[0] |= 2;
				double ef_calc = Straight(
					(he_d1_energy - t0_param[0][0]) / t0_param[0][1],
					(calc_he_d2_energy - t0_param[1][0]) / t0_param[1][1],
					sa12, sb12
				);
				pid_d1d2_calc_straight.Fill(ef_calc);
				if (ef_calc > 6100 && ef_calc < 6600) straight[1] |= 2;
			}
			if (event.layer[0] == 2) {
				double ef = Straight(
					(be_d2_energy - t0_param[1][0]) / t0_param[1][1],
					(be_d3_energy - t0_param[2][0]) / t0_param[2][1],
					sa23, sb23
				);
				pid_d2d3_straight.Fill(ef);
				if (ef > 24000 && ef < 25000) straight[0] |= 1;
				double ef_calc = Straight(
					(calc_be_d2_energy - t0_param[1][0]) / t0_param[1][1],
					(be_d3_energy - t0_param[2][0]) / t0_param[2][1],
					sa23, sb23
				);
				pid_d2d3_calc_straight.Fill(ef_calc);
				if (ef_calc > 24000 && ef_calc < 25000) straight[1] |= 1;
			}
			if (event.layer[1] == 2) {
				double ef = Straight(
					(he_d2_energy - t0_param[1][0]) / t0_param[1][1],
					(he_d3_energy - t0_param[2][0]) / t0_param[2][1],
					sa23, sb23
				);
				pid_d2d3_straight.Fill(ef);
				if (ef > 7000 && ef < 8000) straight[0] |= 2;
				double ef_calc = Straight(
					(calc_he_d2_energy - t0_param[1][0]) / t0_param[1][1],
					(he_d3_energy - t0_param[2][0]) / t0_param[2][1],
					sa23, sb23
				);
				pid_d2d3_calc_straight.Fill(ef_calc);
				if (ef_calc > 7000 && ef_calc < 8000) straight[1] |= 2;
			}
			if (event.layer[1] > 2) {
				straight[0] |= 2;
				straight[1] |= 2;
			}

			// fill Q values
			for (size_t i = 0; i < 4; ++i) hq[i].Fill(q[i]);
		}

		// fill to tree
		opt.Fill();
	}

	// fit D1D2-straight
	TF1 *fd1d2 = new TF1("fd1d2", "gaus(0)+gaus(3)", 5000, 25000);
	constexpr double d1d2_init_param[6] = {
		100, 6500, 100,
		200, 21000, 100
	};
	fd1d2->SetParameters(d1d2_init_param);
	fd1d2->SetParLimits(2, 0, 1000);
	fd1d2->SetParLimits(5, 0, 1000);
	fd1d2->SetNpx(1000);
	pid_d1d2_straight.Fit(fd1d2, "RQ+");
	std::cout << "D1D2 straight "
		<< fd1d2->GetParameter(1) << ", " << fd1d2->GetParameter(2) << ", "
		<< fd1d2->GetParameter(4) << ", " << fd1d2->GetParameter(5) << "\n";

	// fit D1D2-calculated-straight
	TF1 *fd1d2c = new TF1("fd1d2c", "gaus(0)+gaus(3)", 5000, 25000);
	fd1d2c->SetParameters(d1d2_init_param);
	fd1d2c->SetNpx(1000);
	fd1d2c->SetParLimits(2, 0, 1000);
	fd1d2c->SetParLimits(5, 0, 1000);
	pid_d1d2_calc_straight.Fit(fd1d2c, "RQ+");
	std::cout << "D1D2 calculated straight "
		<< fd1d2c->GetParameter(1) << ", " << fd1d2c->GetParameter(2) << ", "
		<< fd1d2c->GetParameter(4) << ", " << fd1d2c->GetParameter(5) << "\n";

	// fit D2D3-straight
	TF1 *fd2d3 = new TF1("fd2d3", "gaus(0)+gaus(3)", 5000, 30000);
	constexpr double d2d3_init_param[6] = {
		100, 7500, 150,
		20, 24500, 300
	};
	fd2d3->SetParameters(d2d3_init_param);
	fd2d3->SetParLimits(2, 0, 1000);
	fd2d3->SetParLimits(5, 0, 1000);
	fd2d3->SetNpx(1000);
	pid_d2d3_straight.Fit(fd2d3, "RQ+");
	std::cout << "D2D3 straight "
		<< fd2d3->GetParameter(1) << ", " << fd2d3->GetParameter(2) << ", "
		<< fd2d3->GetParameter(4) << ", " << fd2d3->GetParameter(5) << "\n";

	// fit D2D3-calculated-straight
	TF1 *fd2d3c = new TF1("fd2d3c", "gaus(0)+gaus(3)", 5000, 30000);
	constexpr double d2d3_calc_init_param[6] = {
		200, 7500, 100,
		60, 24500, 300
	};
	fd2d3c->SetParameters(d2d3_calc_init_param);
	fd2d3c->SetParLimits(2, 0, 1000);
	fd2d3c->SetParLimits(5, 0, 1000);
	fd2d3c->SetNpx(1000);
	pid_d2d3_calc_straight.Fit(fd2d3c, "RQ+");
	std::cout << "D2D3 calculated straight "
		<< fd2d3c->GetParameter(1) << ", " << fd2d3c->GetParameter(2) << ", "
		<< fd2d3c->GetParameter(4) << ", " << fd2d3c->GetParameter(5) << "\n";

	// fit Q0 spectrum
	TF1 *fq0 = new TF1("fq0", "gaus(0)+gaus(3)+gaus(6)", -23, -8);
	constexpr double q_init_param[9] = {
		50, -12, 2,
		180, -15.5, 2,
		140, -18, 2
	};
	fq0->SetParameters(q_init_param);
	fq0->SetParLimits(2, 0, 10);
	fq0->SetParLimits(5, 0, 10);
	fq0->SetParLimits(8, 0, 10);
	hq[0].Fit(fq0, "QR+");
	std::cout << "Q0 "
		<< fq0->GetParameter(1) << ", " << fq0->GetParameter(2) << ", "
		<< fq0->GetParameter(4) << ", " << fq0->GetParameter(5) << ", "
		<< fq0->GetParameter(7) << ", " << fq0->GetParameter(8) << "\n";

	// fit Q2 spectrum
	TF1 *fq2 = new TF1("fq2", "gaus(0)+gaus(3)+gaus(6)", -23, -8);
	fq2->SetParameters(q_init_param);
	fq2->SetParLimits(2, 0, 10);
	fq2->SetParLimits(5, 0, 10);
	fq2->SetParLimits(8, 0, 10);
	hq[2].Fit(fq2, "QR+");
	std::cout << "Q2 "
		<< fq2->GetParameter(1) << ", " << fq2->GetParameter(2) << ", "
		<< fq2->GetParameter(4) << ", " << fq2->GetParameter(5) << ", "
		<< fq2->GetParameter(7) << ", " << fq2->GetParameter(8) << "\n";

	// fit Q3 spectrum
	TF1 *fq3 = new TF1("fq3", "gaus(0)+gaus(3)+gaus(6)", -23, -8);
	fq3->SetParameters(q_init_param);
	fq3->SetParLimits(2, 0, 10);
	fq3->SetParLimits(5, 0, 10);
	fq3->SetParLimits(8, 0, 10);
	hq[3].Fit(fq3, "QR+");
	std::cout << "Q3 "
		<< fq3->GetParameter(1) << ", " << fq3->GetParameter(2) << ", "
		<< fq3->GetParameter(4) << ", " << fq3->GetParameter(5) << ", "
		<< fq3->GetParameter(7) << ", " << fq3->GetParameter(8) << "\n";

	// save
	opf.cd();
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 11; ++j) {
			hist_stateless_excited_energy[i][j].Write();
		}
	}
	// fill spectrums
	pid_d1d2.Write();
	pid_d2d3.Write();
	pid_d1d2_calc.Write();
	pid_d2d3_calc.Write();
	pid_d1d2_straight.Write();
	pid_d2d3_straight.Write();
	pid_d1d2_calc_straight.Write();
	pid_d2d3_calc_straight.Write();
	for (size_t i = 0; i < 4; ++i) hq[i].Write();
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();

	return 0;
}