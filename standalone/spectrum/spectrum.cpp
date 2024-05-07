#include <iostream>
#include <fstream>

#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TString.h>
#include <TTree.h>
#include <Math/Vector3D.h>
#include <THStack.h>

#include "include/event/threebody_info_event.h"
#include "include/ppac_track.h"
#include "include/calculator/d2_energy_calculator.h"

using namespace ribll;


constexpr double tafd_phi_start[6] = {
	117.6, 57.6, -2.4, -62.4, -122.4, 177.6
};

constexpr double theory_q[4] = {
	-12.0125, -15.3805, -18.1915, -19.5545
};


/// @brief rebuild threebody reaction process
/// @param[in] event input event
/// @param[in] csi_energy CsI energy
/// @param[in] tx reaction point x
/// @param[in] ty reaction point y
/// @returns Q value
///
double ThreeBodyProcess(
	const ThreeBodyInfoEvent &event,
	double csi_energy,
	double tafx,
	double tafy
) {
	// target point
	double tx = event.xptx;
	double ty = event.xpty;

	// 10Be momentum
	double be_momentum = MomentumFromKinetic(mass_10be, event.t0_energy[0]);
	// 10Be momentum vector
	ROOT::Math::XYZVector p_be(
		event.be_x[0] - tx,
		event.be_y[0] - ty,
		100.0
	);
	p_be = p_be.Unit() * be_momentum;

	// 4He momentum
	double he_momentum = MomentumFromKinetic(mass_4he, event.t0_energy[1]);
	// 4He momentum vector
	ROOT::Math::XYZVector p_he(
		event.he_x[0] - tx,
		event.he_y[0] - ty,
		100.0
	);
	p_he = p_he.Unit() * he_momentum;

	double taf_energy = event.tafd_energy + csi_energy;
	// 2H momentum
	double d_momentum = MomentumFromKinetic(mass_2h, taf_energy);
	// 2H momentum vector
	ROOT::Math::XYZVector p_d(
		tafx - tx,
		tafy - ty,
		135.0
	);
	p_d = p_d.Unit() * d_momentum;

	// beam 14C momentum vector
	ROOT::Math::XYZVector p_beam = p_be + p_he + p_d;

	// 14C momentum
	double beam_momentum = p_beam.R();
	// 14C kinematic energy
	double c14_kinetic =
		sqrt(pow(beam_momentum, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

	// three-fold Q value
	double q = event.t0_energy[0] + event.t0_energy[1]
		+ taf_energy - c14_kinetic;

	return q;
}


/// @brief calculate 14C exicted energy
/// @param[in] event input event
/// @param[in] state 10Be state, 0-ground state, 1-3.368MeV, 2-6.179MeV, 3-7.542MeV
/// @returns excited energy of 14C
///
double TwoBodyProcess(
	const ThreeBodyInfoEvent &event,
	const int state
) {
	double excited_10be = 0.0;
	if (state == 1) excited_10be = 3.368;
	else if (state == 2) excited_10be = 6.179;
	else if (state == 3) excited_10be = 7.542;

	// 10Be momentum
	double be_momentum = MomentumFromKinetic(
		mass_10be + excited_10be, event.t0_energy[0]
	);
	// 10Be momentum vector
	ROOT::Math::XYZVector p_be(
		event.be_x[0] - event.xptx,
		event.be_y[0] - event.xpty,
		100.0
	);
	p_be = p_be.Unit() * be_momentum;

	// 4He momentum
	double he_momentum = MomentumFromKinetic(mass_4he, event.t0_energy[1]);
	// 4He momentum vector
	ROOT::Math::XYZVector p_he(
		event.he_x[0] - event.xptx,
		event.he_y[0] - event.xpty,
		100.0
	);
	p_he = p_he.Unit() * he_momentum;

	// excited 14C momentum vector
	ROOT::Math::XYZVector p_c = p_be + p_he;
	// excited 14C momentum
	double c_momentum = p_c.R();
	// excited 14C total energy
	double c_energy =
		(event.t0_energy[0] + mass_10be + excited_10be)
		+ (event.t0_energy[1] + mass_4he);
	// excited 14C mass
	double excited_c_mass = sqrt(
		pow(c_energy, 2.0) - pow(c_momentum, 2.0)
	);
	// excited energy of 14C
	double excited_14c = excited_c_mass - mass_14c;

	return excited_14c;
}


/// @brief check the Q values and get the state and type  
/// @param[in] q Q values
/// @param[out] type -1: invalid, 0:three same values,
/// 	1:two same values, 2:different values
/// @returns confirmed state
///
int CheckRepresentQ(double *q, int &type) {
	// check linear Q values
	double diff[3];
	int state[3];
	int counts[4] = {0, 0, 0, 0};
	// search for closest state
	for (int i = 0; i < 3; ++i) {
		diff[i] = q[i] - theory_q[0];
		state[i] = 0;
		for (int j = 1; j < 4; ++j) {
			double t = q[i] - theory_q[j];
			if (fabs(t) < fabs(diff[i])) {
				diff[i] = t;
				state[i] = j;
			}
		}
		// increase counts
		if (q[i] < -11.0 && q[i] > -20.5) {
			++counts[state[i]];
		}
	}
	// maximum counts
	int max_count = counts[0];
	// state with maximum counts
	int max_state = 0;
	for (int i = 1; i < 4; ++i) {
		if (counts[i] > max_count) {
			max_count = counts[i];
			max_state = i;
		}
	}

	type = -1;
	if (max_count == 1) {
		type = 2;
	} else if (max_count == 2) {
		type = 1;
	} else if (max_count == 3) {
		type = 0;
	}

	return max_state;
}


class Spectrum {
public:

	Spectrum(int n): n_(n) {}

	double operator()(double *x, double *par) {
		double result = 0.0;
		for (int i = 0; i < n_; ++i) {
			result += par[3*i] * exp(-0.5 * pow((x[0] - par[3*i+1]) / par[3*i+2], 2.0));
		}
		return result;
	}

private:
	int n_;
};


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

	// simulated file to get efficiency
	TString efficiency_file_name = TString::Format(
		"%s%sefficiency-0002.root",
		kGenerateDataPath,
		kSimulateDir
	);
	// simulate file
	TFile efficiency_file(efficiency_file_name, "read");
	TGraph *geff[3];
	for (int i = 0; i < 3; ++i) {
		geff[i] = (TGraph*)efficiency_file.Get(TString::Format("g%d", i));
	}
	if (!geff[0] || !geff[1] || !geff[2]) {
		std::cerr << "Error: Get efficiency graph from "
			<< efficiency_file_name << " failed.\n";
	}

	// output file name
	TString output_file_name = TString::Format(
		"%s%sthreebody%s.root",
		kGenerateDataPath,
		kSpectrumDir,
		suffix.empty() ? "" : ("-"+suffix).c_str()
	);
	// output file
	TFile output_file(output_file_name, "recreate");

	// spectrum of 14C decay to 10Be ground state
	TH1F c_spec_0(
		"hsc0", "spectrum of 14C to 10Be ground state",
		100, 12, 32
	);
	// spectrum of 14C decay to 10Be 3.368 MeV state
	TH1F c_spec_1(
		"hsc1", "spectrum of 14C to 10Be 3.368MeV state",
		100, 12, 32
	);
	// spectrum of 14C decay to 10Be 6.179 MeV state
	TH1F c_spec_2(
		"hsc2", "spectrum of 14C to 10Be 6.179MeV state",
		100, 12, 32
	);
	// spectrum of 14C decay to 10Be 7.542 MeV state
	TH1F c_spec_3(
		"hsc3", "spectrum of 14C to 10Be 7.542MeV state",
		100, 12, 32
	);
	// spectrum in the range of hjx's article
	TH1F hjx_spectrum0("hjx0", "spectrum", 35, 12, 19);
	// sepectrum of first excited state in hjx's article range
	TH1F hjx_spectrum1("hjx1", "spectrum", 45, 16, 25);
	// spectrum of 6MeV excited state in Baba's article predicted range
	TH1F sigma_predicted_spectrum2("baba2", "spectrum", 55, 19, 30);
	// hjx_spectrum0.SetLineColor(kBlack);
	// hjx_spectrum1.SetLineColor(kBlack);
	// sigma_predicted_spectrum2.SetLineColor(kBlack);

	// spectrum of 14C decay to 10Be ground state, represent events
	TH1F c_rep_spec_0(
		"hrsc0", "spectrum of 14C to 10Be ground state (represent)",
		100, 12, 32
	);
	// spectrum of 14C decay to 10Be 3.368 MeV state, represent events
	TH1F c_rep_spec_1(
		"hrsc1", "spectrum of 14C to 10Be 3.368MeV state (represent)",
		100, 12, 32
	);
	// spectrum of 14C decay to 10Be 6.179 MeV state, represent events
	TH1F c_rep_spec_2(
		"hrsc2", "spectrum of 14C to 10Be 6.179MeV state (represent)",
		100, 12, 32
	);
	// spectrum of 14C decay to 10Be 7.542 MeV state, represent events
	TH1F c_rep_spec_3(
		"hrsc3", "spectrum of 14C to 10Be 7.542MeV state (represent)",
		100, 12, 32
	);
	// spectrum in the range of hjx's article, represent events
	TH1F hjx_rep_spectrum0("hjxr0", "spectrum", 35, 12, 19);
	// sepectrum of 3.368 MeV state in hjx's article range, represent events
	TH1F hjx_rep_spectrum1("hjxr1", "spectrum", 45, 16, 25);
	// spectrum of 6.179 MeV state in Baba's article predicted range, represent
	TH1F sigma_predicted_rep_spectrum2("babar2", "spectrum", 55, 19, 30);
	// hjx_rep_spectrum0.SetLineColor(kBlack);
	// hjx_rep_spectrum1.SetLineColor(kBlack);
	// sigma_predicted_rep_spectrum2.SetLineColor(kBlack);
	// spectrum of 14C decay to 10Be ground state, efficiency corrected
	TH1F c_rep_spec_0_corr(
		"hrsc0c",
		"spectrum of 14C to 10Be ground state (represent, efficiency)",
		100, 12, 32
	);
	// spectrum of 14C decay to 10Be 3.368 MeV state, efficiency corrected
	TH1F c_rep_spec_1_corr(
		"hrsc1c",
		"spectrum of 14C to 10Be 3.368MeV state (represent, efficiency)",
		100, 12, 32
	);
	// spectrum of 14C decay to 10Be 6.179 MeV state, efficiency correted
	TH1F c_rep_spec_2_corr(
		"hrsc2c",
		"spectrum of 14C to 10Be 6.179MeV state (represent, efficiency)",
		100, 12, 32
	);
	// multiplied efficiency to show in the same graph with spectrum  
	TGraph cgeff[3];
	


	// output tree
	TTree opt("tree", "spectrum");
	// output data
	double linear_q[3], power_q[3], opt_q[3];
	double calc_be_q[3], calc_he_q[3], calc_behe_q[3];
	// traditional state
	int linear_state, power_state, opt_state;
	// represent state
	int linear_represent_state, power_represent_state, opt_represent_state;
	int calc_be_represent_state, calc_he_represent_state, calc_behe_represent_state;
	int linear_type, power_type, opt_type;
	int calc_be_type, calc_he_type, calc_behe_type;
	int be_state;
	double c_excited, c_possible_excited[4];
	int ppac_flag, taf_flag, bind, hole;
	// setup output branches
	opt.Branch("linear_q", linear_q, "lq[3]/D");
	opt.Branch("power_q", power_q, "pq[3]/D");
	opt.Branch("opt_q", opt_q, "oq[3]/D");
	opt.Branch("cbe_q", calc_be_q, "c1q[3]/D");
	opt.Branch("che_q", calc_he_q, "c2q[3]/D");
	opt.Branch("cbehe_q", calc_behe_q, "c3q[3]/D");
	opt.Branch("linear_rep_state", &linear_represent_state, "lrstate/I");
	opt.Branch("power_rep_state", &power_represent_state, "prstate/I");
	opt.Branch("opt_rep_state", &opt_represent_state, "orstate/I");
	opt.Branch("cbe_rep_state", &calc_be_represent_state, "c1state/I");
	opt.Branch("che_rep_state", &calc_he_represent_state, "c2state/I");
	opt.Branch("cbehe_rep_state", &calc_behe_represent_state, "c3state/I");
	opt.Branch("linear_type", &linear_type, "ltype/I");
	opt.Branch("power_type", &power_type, "ptype/I");
	opt.Branch("opt_type", &opt_type, "otype/I");
	opt.Branch("cbe_type", &calc_be_type, "c1type/I");
	opt.Branch("che_type", &calc_he_type, "c2type/I");
	opt.Branch("cbehe_type", &calc_behe_type, "c3type/I");
	opt.Branch("linear_state", &linear_state, "lstate/I");	
	opt.Branch("power_state", &power_state, "pstate/I");
	opt.Branch("opt_state", &opt_state, "ostate/I");
	opt.Branch("be_state",  &be_state, "state/I");
	opt.Branch("c_excited", &c_excited, "cex/D");
	opt.Branch("c_possible_excited", c_possible_excited, "cpex[4]/D");
	opt.Branch("ppac_flag", &ppac_flag, "pflag/I");
	opt.Branch("taf_flag", &taf_flag, "tflag/I");
	opt.Branch("bind", &bind, "bind/I");
	opt.Branch("hole", &hole, "hole/I");


	// Q value correct
	constexpr double q_correct[12] = {
		-0.46, -0.25, 0.18, 0.49,
		0.49, 0.55, 0.28, 0.30,
		-0.04, -0.10, -0.12, -0.64
	};

	// calibrate parameters
	double linear_param[12][2];
	double power_param[12][3];

	// read linear calibrate parameters
	// loop TAFs
	for (int index = 0; index < 6; ++index) {
		// file name   
		TString file_name = TString::Format(
			"%s%staf%dcsi-cali-param-2H-linear.txt",
			kGenerateDataPath,
			kCalibrationDir,
			index
		);
		// input stream
		std::ifstream linear_fin(file_name.Data());
		if (!linear_fin.good()) {
			std::cerr << "Error: Get parameters from file "
				<< file_name << " failed.\n";
			return -1;
		}
		linear_fin >> linear_param[index*2][0] >> linear_param[index*2][1]
			>> linear_param[index*2+1][0] >> linear_param[index*2+1][1];
		// close files
		linear_fin.close();
	}
	// show read parameteres
	std::cout << "Read linear calibrate parameters\n";
	for (int i = 0; i < 12; ++i) {
		std::cout << linear_param[i][0] << ", "
			<< linear_param[i][1] << "\n";
	}

	// read power calibrate parameters
	// loop TAFs
	for (int index = 0; index < 6; ++index) {
		// file name   
		TString file_name = TString::Format(
			"%s%staf%dcsi-cali-param-2H.txt",
			kGenerateDataPath,
			kCalibrationDir,
			index
		);
		// input stream
		std::ifstream power_fin(file_name.Data());
		if (!power_fin.good()) {
			std::cerr << "Error: Get parameters from file "
				<< file_name << " failed.\n";
			return -1;
		}
		for (int i = 0; i < 2; ++i) {
			power_fin >> power_param[index*2+i][0]
				>> power_param[index*2+i][1]
				>> power_param[index*2+i][2];
		}
		// close files
		power_fin.close();
	}
	// show read parameteres
	std::cout << "Read power calibrate parameters 2\n";
	for (int i = 0; i < 12; ++i) {
		std::cout << power_param[i][0] << ", "
			<< power_param[i][1] << ", "
			<< power_param[i][2] << "\n";
	}

	elc::D2EnergyCalculator be_calculator("10Be");
	elc::D2EnergyCalculator he_calculator("4He");


	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		// get status
		taf_flag = event.taf_flag;
		ppac_flag = event.ppac_flag;
		bind = event.bind;
		hole = 0;
		hole |= event.hole[0] ? 0x1 : 0;
		hole |= event.hole[1] ? 0x2 : 0;

		// backup kinetic energy
		const double be_energy = event.t0_energy[0];
		const double he_energy = event.t0_energy[1];

		// calculate 10Be energy
		// energy correct
		double be_correct = 0.0;
		double be_d2_energy =
			t0_param[1][0] + t0_param[1][1] * event.be_channel[1];
		if (event.layer[0] == 1) {
			double d1_energy =
				t0_param[0][0] + t0_param[0][1] * event.be_channel[0];
			// stop in T0D2
			double calc_be_energy =
				be_calculator.Energy(0, d1_energy, 0.0, true);
			// get energy correct
			be_correct = calc_be_energy - be_d2_energy;
		} else if (event.layer[0] == 2) {
			double d3_energy =
				t0_param[2][0] + t0_param[2][1] * event.be_channel[2];
			// stop in T0D3
			double calc_be_energy =
				be_calculator.DeltaEnergy(1, d3_energy, 0.0, true);
			// get energy correct
			be_correct = calc_be_energy - be_d2_energy;
		}
		// calculate 4He energy
		double he_correct = 0.0;
		double he_d2_energy =
			t0_param[1][0] + t0_param[1][1] * event.he_channel[1];
		if (event.layer[1] == 1) {
			double d1_energy =
				t0_param[0][0] + t0_param[0][1] * event.he_channel[0];
			// stop in T0D2
			double calc_he_energy =
				he_calculator.Energy(0, d1_energy, 0.0, true);
			// get energy correct
			he_correct = calc_he_energy - he_d2_energy;
		} else if (event.layer[1] >= 2) {
			double d3_energy =
				t0_param[2][0] + t0_param[2][1] * event.he_channel[2];
			// stop in or pass T0D3
			double calc_he_energy = he_calculator.DeltaEnergy(
				1, d3_energy, 0.0, event.layer[1]==2
			);
			// get energy correct
			he_correct = calc_he_energy - he_d2_energy;
		}


		// loop to get different TAF X and Y
		for (int i = 0; i < 3; ++i) {
			// change 2H position
			double ctaf_r = 102.5/16.0 * (event.d_x_strip + i*0.5) + 32.6;
			double ctaf_phi =
				-55.2/8.0 * (event.d_y_strip+0.5)
				+ tafd_phi_start[event.csi_index/2];
			ctaf_phi *= TMath::DegToRad();
			double mid_phi =
				(tafd_phi_start[event.csi_index/2] - 27.6) * TMath::DegToRad();
			// get changed position
			double ctafx = ctaf_r * cos(ctaf_phi) + 34.4*cos(mid_phi);
			double ctafy = ctaf_r * sin(ctaf_phi) + 34.4*sin(mid_phi);

			// get Q values
			// linear CsI energy
			double linear_csi_energy =
				linear_param[event.csi_index][0]
				+ linear_param[event.csi_index][1] * event.csi_channel;
			// linear Q value
			linear_q[i] = ThreeBodyProcess(event, linear_csi_energy, ctafx, ctafy);
			linear_q[i] -= q_correct[event.csi_index];

			// power CsI energy
			double power_csi_energy = pow(
				(event.csi_channel - power_param[event.csi_index][2])
					/ power_param[event.csi_index][0],
				1.0 / power_param[event.csi_index][1]
			);
			// power Q value
			power_q[i] = ThreeBodyProcess(event, power_csi_energy, ctafx, ctafy);
			power_q[i] -= q_correct[event.csi_index];

			// optimized CsI energy
			double opt_csi_energy = pow(
				(event.csi_channel - csi_param[event.csi_index][2])
					/ csi_param[event.csi_index][0],
				1.0 / csi_param[event.csi_index][1] 
			);
			// optimized Q value
			opt_q[i] = ThreeBodyProcess(event, opt_csi_energy, ctafx, ctafy);
			opt_q[i] -= q_correct[event.csi_index];


			// calculate Be energy
			event.t0_energy[0] = be_energy + be_correct;
			calc_be_q[i] = ThreeBodyProcess(event, power_csi_energy, ctafx, ctafy);
			calc_be_q[i] -= q_correct[event.csi_index];

			// calculate BeHe energy
			event.t0_energy[1] = he_energy + he_correct;
			calc_behe_q[i] = ThreeBodyProcess(event, power_csi_energy, ctafx, ctafy);
			calc_behe_q[i] -= q_correct[event.csi_index];

			// calculate He energy
			event.t0_energy[0] = be_energy;
			calc_he_q[i] = ThreeBodyProcess(event, power_csi_energy, ctafx, ctafy);
			calc_he_q[i] -= q_correct[event.csi_index];

			// backup energy
			event.t0_energy[1] = he_energy;
		}


		// get state and type
		linear_represent_state = CheckRepresentQ(linear_q, linear_type);
		power_represent_state = CheckRepresentQ(power_q, power_type);
		opt_represent_state = CheckRepresentQ(opt_q, opt_type);
		calc_be_represent_state = CheckRepresentQ(calc_be_q, calc_be_type);
		calc_he_represent_state = CheckRepresentQ(calc_he_q, calc_he_type);
		calc_behe_represent_state = CheckRepresentQ(calc_behe_q, calc_behe_type);


		if (linear_q[1] < -11 && linear_q[1] > -13.5) linear_state = 0;
		else if (linear_q[1] < -14.5 && linear_q[1] > -16.5) linear_state = 1;
		else if (linear_q[1] < -17 && linear_q[1] > -19) linear_state = 2;
		else if (linear_q[1] < -19 && linear_q[1] > -20.5) linear_state = 3;
		else linear_state = -1;

		if (power_q[1] < -11.5 && power_q[1] > -12.5) power_state = 0;
		else if (power_q[1] < -15 && power_q[1] > -16.2) power_state = 1;
		else if (power_q[1] < -17.5 && power_q[1] > -18.5) power_state = 2;
		else if (power_q[1] < -19 && power_q[1] > -20.5) power_state = 3;
		else power_state = -1;

		if (opt_q[1] < -11 && opt_q[1] > -13.5) opt_state = 0;
		else if (opt_q[1] < -14.5 && opt_q[1] > -16.5) opt_state = 1;
		else if (opt_q[1] < -17 && opt_q[1] > -19) opt_state = 2;
		else if (opt_q[1] < -19 && opt_q[1] > -20.5) opt_state = 3;
		else opt_state = -1;

		be_state = power_represent_state;
		c_excited = TwoBodyProcess(event, be_state);
		
		for (int i = 0; i < 4; ++i) {
			c_possible_excited[i] = TwoBodyProcess(event, i);
		}

		if (
			taf_flag == 0
			&& (ppac_flag & 1) == 1
			&& bind == 0
			&& hole == 0
		) {
			if (power_state == 0){
				c_spec_0.Fill(c_excited);
				hjx_spectrum0.Fill(c_excited);
			} else if (power_state == 1) {
				c_spec_1.Fill(c_excited);
				hjx_spectrum1.Fill(c_excited);
			} else if (power_state == 2) {
				c_spec_2.Fill(c_excited);
				sigma_predicted_spectrum2.Fill(c_excited);
			} else if (power_state == 3) {
				c_spec_3.Fill(c_excited);
			}

			if (be_state == 0){
				c_rep_spec_0.Fill(c_excited);
				hjx_rep_spectrum0.Fill(c_excited);
			} else if (be_state == 1) {
				c_rep_spec_1.Fill(c_excited);
				hjx_rep_spectrum1.Fill(c_excited);
			} else if (be_state == 2) {
				c_rep_spec_2.Fill(c_excited);
				sigma_predicted_rep_spectrum2.Fill(c_excited);
			} else if (be_state == 3) {
				c_rep_spec_3.Fill(c_excited);
			}
		}

		// fill to tree
		opt.Fill();
	}

	// efficiency corrected spectrum
	TH1F *spec_rep[3] = {&c_spec_0, &c_spec_1, &c_spec_2};
	TH1F *spec_corr[3] = {
		&c_rep_spec_0_corr, &c_rep_spec_1_corr, &c_rep_spec_2_corr
	};
	for (int i = 0; i < 3; ++i) {
		int bins = spec_rep[i]->GetNbinsX();
		for (int bin = 1; bin <= bins; ++bin) {
			spec_corr[i]->SetBinContent(
				bin,
				spec_rep[i]->GetBinContent(bin)
				/ geff[i]->Eval(spec_rep[i]->GetBinCenter(bin)) 
			);
		}
	}


	// fit hjx range spectrum
	Spectrum hjx_peaks(6);
	// fitting function
	TF1 fit_hjx("fhjx", hjx_peaks, 12, 19, 21);
	double fit_hjx_initial_parameters[18] = {
		10.0, 13.6, 0.1,
		20.0, 14.9, 0.05,
		20.0, 15.6, 0.05,
		10.0, 16.4, 0.1,
		10.0, 17.3, 0.1,
		20.0, 18.2, 0.4
	};

	fit_hjx.SetParameters(fit_hjx_initial_parameters);

	// fit_hjx.SetParLimits(0, 1.0, 3.0);
	fit_hjx.SetParLimits(1, 13.4, 14.0);
	fit_hjx.SetParLimits(2, 0.0, 1.0);

	// fit_hjx.SetParLimits(3, 2.0, 10.0);
	fit_hjx.SetParLimits(4, 14.5, 15.2);
	fit_hjx.SetParLimits(5, 0.0, 0.4);

	// fit_hjx.SetParLimits(6, 2.0, 30.0);
	fit_hjx.SetParLimits(7, 15.2, 16.0);
	fit_hjx.SetParLimits(8, 0.0, 0.2);

	// fit_hjx.SetParLimits(9, 0.0, 10.0);
	fit_hjx.SetParLimits(10, 15.9, 16.7);
	fit_hjx.SetParLimits(11, 0.0, 0.25);

	// fit_hjx.SetParLimits(12, 0.0, 10.0);
	fit_hjx.SetParLimits(13, 16.7, 17.6);
	fit_hjx.SetParLimits(14, 0.0, 0.15);

	// fit_hjx.SetParLimits(15, 0.0, 10.0);
	fit_hjx.SetParLimits(16, 17.6, 18.6);
	fit_hjx.SetParLimits(17, 0.0, 1.0);

	fit_hjx.SetNpx(1000);
	hjx_rep_spectrum0.Fit(&fit_hjx, "R+");
	double hjx_final_parameters[18];
	fit_hjx.GetParameters(hjx_final_parameters);
	TF1 *fit_hjxs[6];
	for (int i = 0; i < 6; ++i) {
		fit_hjxs[i] = new TF1(
			TString::Format("fhjx0%d", i), "gaus", 12, 19
		);
		fit_hjxs[i]->SetParameters(hjx_final_parameters+3*i);
		fit_hjxs[i]->SetLineColor(kBlue);
		fit_hjxs[i]->SetNpx(200);
		hjx_rep_spectrum0.GetListOfFunctions()->Add(fit_hjxs[i]);
	}

	std::cout << "Fit hjx spectrum result:\n";
	for (int i = 0; i < 6; ++i) {
		std::cout << hjx_final_parameters[i*3] << " " << hjx_final_parameters[i*3+1]
			<< " " << hjx_final_parameters[i*3+2] << "\n";
	}



	// // first excited state
	// Spectrum first_peaks(3);
	// // fitting function
	// TF1 fit_hjx1("fhjx1", first_peaks, 16, 19, 9);
	// double fit_hjx1_initial_parameters[9] = {
	// 	10.0, 17.3, 0.08,
	// 	5.0, 17.9, 0.1,
	// 	15.0, 18.5, 0.1,
	// };

	// fit_hjx1.SetParameters(fit_hjx1_initial_parameters);

	// fit_hjx1.SetParLimits(0, 2.0, 20.0);
	// fit_hjx1.SetParLimits(1, 17.0, 17.5);
	// fit_hjx1.SetParLimits(2, 0.0, 0.2);

	// fit_hjx1.SetParLimits(3, 2.0, 20.0);
	// fit_hjx1.SetParLimits(4, 17.2, 18.0);
	// fit_hjx1.SetParLimits(5, 0.0, 0.4);

	// fit_hjx1.SetParLimits(6, 2.0, 30.0);
	// fit_hjx1.SetParLimits(7, 18.0, 19.0);
	// fit_hjx1.SetParLimits(8, 0.0, 0.4);


	// fit_hjx1.SetNpx(1000);
	// hjx_spectrum1.Fit(&fit_hjx1, "R+");
	// double hjx1_final_parameters[9];
	// fit_hjx1.GetParameters(hjx1_final_parameters);
	// TF1 *fit_hjxs1[3];
	// for (int i = 0; i < 3; ++i) {
	// 	fit_hjxs1[i] = new TF1(
	// 		TString::Format("fhjx1%d", i), "gaus", 16, 19
	// 	);
	// 	fit_hjxs1[i]->SetParameters(hjx1_final_parameters+3*i);
	// 	fit_hjxs1[i]->SetLineColor(kBlue);
	// 	fit_hjxs1[i]->SetNpx(200);
	// 	hjx_spectrum1.GetListOfFunctions()->Add(fit_hjxs1[i]);
	// }

	// std::cout << "Fit hjx spectrum first excited result:\n";
	// for (int i = 0; i < 3; ++i) {
	// 	std::cout << hjx1_final_parameters[i*3] << " " << hjx1_final_parameters[i*3+1]
	// 		<< " " << hjx1_final_parameters[i*3+2] << "\n";
	// }



	// // 6MeV excited state
	// Spectrum second_peaks(3);
	// // fitting function
	// TF1 fit_baba2("fbaba2", second_peaks, 21, 30, 9);
	// double fit_baba2_initial_parameters[9] = {
	// 	20.0, 21.4, 0.05,
	// 	20.0, 22.2, 0.05,
	// 	10.0, 23.6, 0.2
	// };

	// fit_baba2.SetParameters(fit_baba2_initial_parameters);

	// fit_baba2.SetParLimits(0, 2.0, 50.0);
	// fit_baba2.SetParLimits(1, 21.0, 22.0);
	// fit_baba2.SetParLimits(2, 0.0, 1.0);

	// fit_baba2.SetParLimits(3, 2.0, 50.0);
	// fit_baba2.SetParLimits(4, 21.8, 23.0);
	// fit_baba2.SetParLimits(5, 0.0, 1.0);

	// fit_baba2.SetParLimits(6, 0.0, 50.0);
	// fit_baba2.SetParLimits(7, 23.0, 25.0);
	// fit_baba2.SetParLimits(8, 0.0, 1.0);


	// fit_baba2.SetNpx(1000);
	// sigma_predicted_spectrum2.Fit(&fit_baba2, "R+");
	// double baba2_final_parameters[9];
	// fit_baba2.GetParameters(baba2_final_parameters);
	// TF1 *fit_babas2[3];
	// for (int i = 0; i < 3; ++i) {
	// 	fit_babas2[i] = new TF1(
	// 		TString::Format("fbaba2%d", i), "gaus", 21, 30
	// 	);
	// 	fit_babas2[i]->SetParameters(baba2_final_parameters+3*i);
	// 	fit_babas2[i]->SetLineColor(kBlue);
	// 	fit_babas2[i]->SetNpx(200);
	// 	sigma_predicted_spectrum2.GetListOfFunctions()->Add(fit_babas2[i]);
	// }

	// std::cout << "Fit predicted sigma bond spectrum result:\n";
	// for (int i = 0; i < 3; ++i) {
	// 	std::cout << baba2_final_parameters[i*3] << " " << baba2_final_parameters[i*3+1]
	// 		<< " " << baba2_final_parameters[i*3+2] << "\n";
	// }

	const double energy_base[3] = {12.02, 15.39, 18.20};
	const double ratio[3] = {100.0, 150.0, 150.0};
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 100; ++j) {
			double e = energy_base[i] + 0.2 * j;
			cgeff[i].AddPoint(
				e, geff[i]->Eval(e) * ratio[i] * 0.5
			);
		}
	}

	THStack c_rep_stack("hrs", "spectrum of 14C (represent)");
	c_rep_stack.Add(&c_rep_spec_0);
	c_rep_stack.Add(&c_rep_spec_1);
	c_rep_stack.Add(&c_rep_spec_2);

	// save
	output_file.cd();
	c_spec_0.Write();
	c_spec_1.Write();
	c_spec_2.Write();
	c_spec_3.Write();
	hjx_spectrum0.Write();
	hjx_spectrum1.Write();
	sigma_predicted_spectrum2.Write();
	c_rep_spec_0.Write();
	c_rep_spec_1.Write();
	c_rep_spec_2.Write();
	c_rep_spec_3.Write();
	c_rep_spec_0.SetLineColor(kBlue);
	c_rep_spec_1.SetLineColor(kBlack);
	c_rep_spec_2.SetLineColor(kRed);
	c_rep_stack.Write();
	c_rep_spec_0_corr.Write();
	c_rep_spec_1_corr.Write();
	c_rep_spec_2_corr.Write();
	hjx_rep_spectrum0.Write();
	hjx_rep_spectrum1.Write();
	sigma_predicted_rep_spectrum2.Write();
	for (int i = 0; i < 3; ++i) {
		cgeff[i].Write(TString::Format("gef%d", i));
	}
	// save tree
	opt.Write();
	// close files
	output_file.Close();
	return 0;
}