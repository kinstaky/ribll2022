#include <cmath>
#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TH1F.h>
#include <TString.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include <ceres/ceres.h>
#include <glog/logging.h>

#include "include/event/threebody_info_event.h"
#include "include/calculator/delta_energy_calculator.h"
#include "include/ppac_track.h"

// #define THRES_OPT
#define POS_OPT
// #define PPAC_OPT

using namespace ribll;


double ThreeBodyProcess(
	const ThreeBodyInfoEvent &event,
	const double *csi_opt_param,
	const double *pos_opt_param,
	double &d_kinetic,
	double &c_kinetic
) {
	// // 10Be kinetic energy
	// // T0D1 energy
	// be_kinetic = t0_param[0][0] + t0_param[0][1] * event.be_channel[0];
	// // T0D2 energy
	// if (calculated_t0d2) {
	// 	be_kinetic += be_calc->Energy(0, be_kinetic);
	// } else {
	// 	be_kinetic += t0_param[1][0] + t0_param[1][1] * event.be_channel[1];
	// }
	// // T0D3 energy
	// if (event.layer[0] > 1) {
	// 	be_kinetic += t0_param[2][0] + t0_param[2][1] * event.be_channel[2];
	// }
	// // T0S1 energy
	// if (event.layer[0] > 2) {
	// 	be_kinetic += t0_param[3][0] + t0_param[3][1] * event.ssd_channel[0];
	// }
	// // T0S2 energy
	// if (event.layer[0] > 3) {
	// 	be_kinetic += t0_param[4][0] + t0_param[4][1] * event.ssd_channel[1];
	// }
	// // T0S3 energy
	// if (event.layer[0] > 4) {
	// 	be_kinetic += t0_param[5][0] + t0_param[5][1] * event.ssd_channel[2];
	// }

	// // 4He kinetic energy
	// // T0D1 energy
	// he_kinetic = t0_param[0][0] + t0_param[0][1] * event.he_channel[0];
	// // T0D2 energy
	// if (calculated_t0d2) {
	// 	he_kinetic += he_calc->Energy(0, he_kinetic);
	// } else {
	// 	he_kinetic += t0_param[1][0] + t0_param[1][1] * event.he_channel[1];
	// }
	// // T0D3 energy
	// if (event.layer[1] > 1) {
	// 	he_kinetic += t0_param[2][0] + t0_param[2][1] * event.he_channel[2];
	// }
	// // T0S1 energy
	// if (event.layer[1] > 2) {
	// 	he_kinetic += t0_param[3][0] + t0_param[3][1] * event.ssd_channel[0];
	// }
	// // T0S2 energy
	// if (event.layer[1] > 3) {
	// 	he_kinetic += t0_param[4][0] + t0_param[4][1] * event.ssd_channel[1];
	// }
	// // T0S3 energy
	// if (event.layer[1] > 4) {
	// 	he_kinetic += t0_param[5][0] + t0_param[5][1] * event.ssd_channel[2];
	// }

	// 2H kinetic energy
	d_kinetic = event.tafd_energy + pow(
		(event.csi_channel - csi_opt_param[2]) / csi_opt_param[0],
		1.0 / csi_opt_param[1]
	);


	// double ppac_correct[2][3] = {
	// 	{ppac_param[0], 0.0, 0.0},
	// 	{ppac_param[1], 0.0, 0.0}
	// };
	// PpacOffsetX(ppac_correct[0][0], ppac_correct[0][1], ppac_correct[0][2]);
	// PpacOffsetY(ppac_correct[1][0], ppac_correct[1][1], ppac_correct[1][2]);
	// // calculate reaction point
	// double ppac_cx[3], ppac_cy[3];
	// for (int i = 0; i < 3; ++i) {
	// 	ppac_cx[i] = event.ppac_x[i] - ppac_correct[0][i];
	// 	ppac_cy[i] = event.ppac_y[i] - ppac_correct[1][i];
	// }
	// // slope and intercept
	// double xk, yk, xb, yb;
	// TrackMultiplePpac(event.ppac_xflag, ppac_xz, ppac_cx, xk, xb);
	// TrackMultiplePpac(event.ppac_yflag, ppac_yz, ppac_cy, yk, yb);

	// 10Be momentum
	double be_momentum = MomentumFromKinetic(mass_10be, event.t0_energy[0]);
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

	// 2H momentum
	double d_momentum = MomentumFromKinetic(mass_2h, d_kinetic);
	// 2H momentum vector
	ROOT::Math::XYZVector p_d(
		event.d_x - event.xptx + pos_opt_param[0],
		event.d_y - event.xpty + pos_opt_param[1],
		135.0
	);
	p_d = p_d.Unit() * d_momentum;

	// beam 14C momentum vector
	ROOT::Math::XYZVector p_c = p_be + p_he + p_d;
// std::cout << p_be.X() << ", " << p_be.Y() << ", " << p_be.Z() << "\n"
// 	<< p_he.X() << ", " << p_he.Y() << ", " << p_he.Z() << "\n"
// 	<< p_d.X() << ", " << p_d.Y() << ", " << p_d.Z() << "\n"
// 	<< p_c.X() << ", " << p_c.Y() << ", " << p_c.Z() << "\n";

	// 14C momentum
	double c_momentum = p_c.R();
	// 14C kinetic energy
	c_kinetic =
		sqrt(pow(c_momentum, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

	// three-fold Q value
	double q = event.t0_energy[0] + event.t0_energy[1] + d_kinetic - c_kinetic;

	return q;
}


// double TwoBodyProcess(
// 	const ThreeBodyInfoEvent &event,
// 	const double *ppac_param,
// 	double excited_10be
// ) {
// 	// 10Be kinetic energy
// 	// T0D1 energy
// 	double be_kinetic = t0_param[0][0] + t0_param[0][1] * event.be_channel[0];
// 	// T0D2 energy
// 	be_kinetic += t0_param[1][0] + t0_param[1][1] * event.be_channel[1];
// 	// T0D3 energy
// 	if (event.layer[0] > 1) {
// 		be_kinetic += t0_param[2][0] + t0_param[2][1] * event.be_channel[2];
// 	}
// 	// T0S1 energy
// 	if (event.layer[0] > 2) {
// 		be_kinetic += t0_param[3][0] + t0_param[3][1] * event.ssd_channel[0];
// 	}
// 	// T0S2 energy
// 	if (event.layer[0] > 3) {
// 		be_kinetic += t0_param[4][0] + t0_param[4][1] * event.ssd_channel[1];
// 	}
// 	// T0S3 energy
// 	if (event.layer[0] > 4) {
// 		be_kinetic += t0_param[5][0] + t0_param[5][1] * event.ssd_channel[2];
// 	}

// 	// 4He kinetic energy
// 	// T0D1 energy
// 	double he_kinetic = t0_param[0][0] + t0_param[0][1] * event.he_channel[0];
// 	// T0D2 energy
// 	he_kinetic += t0_param[1][0] + t0_param[1][1] * event.he_channel[1];
// 	// T0D3 energy
// 	if (event.layer[1] > 1) {
// 		he_kinetic += t0_param[2][0] + t0_param[2][1] * event.he_channel[2];
// 	}
// 	// T0S1 energy
// 	if (event.layer[1] > 2) {
// 		he_kinetic += t0_param[3][0] + t0_param[3][1] * event.ssd_channel[0];
// 	}
// 	// T0S2 energy
// 	if (event.layer[1] > 3) {
// 		he_kinetic += t0_param[4][0] + t0_param[4][1] * event.ssd_channel[1];
// 	}
// 	// T0S3 energy
// 	if (event.layer[1] > 4) {
// 		he_kinetic += t0_param[5][0] + t0_param[5][1] * event.ssd_channel[2];
// 	}

// 	double ppac_correct[2][3] = {
// 		{ppac_param[0], 0.0, 0.0},
// 		{ppac_param[1], 0.0, 0.0}
// 	};
// 	PpacOffsetX(ppac_correct[0][0], ppac_correct[0][1], ppac_correct[0][2]);
// 	PpacOffsetY(ppac_correct[1][0], ppac_correct[1][1], ppac_correct[1][2]);
// 	// calculate reaction point
// 	double ppac_cx[3], ppac_cy[3];
// 	for (int i = 0; i < 3; ++i) {
// 		ppac_cx[i] = event.ppac_x[i] - ppac_correct[0][i];
// 		ppac_cy[i] = event.ppac_y[i] - ppac_correct[1][i];
// 	}
// 	// slope and intercept
// 	double xk, yk, xb, yb;
// 	TrackMultiplePpac(event.ppac_xflag, ppac_xz, ppac_cx, xk, xb);
// 	TrackMultiplePpac(event.ppac_yflag, ppac_yz, ppac_cy, yk, yb);

// 	// 10Be momentum
// 	double be_momentum = MomentumFromKinetic(
// 		mass_10be + excited_10be, be_kinetic
// 	);
// 	// 10Be momentum vector
// 	ROOT::Math::XYZVector p_be(
// 		event.be_x[0] - xb,
// 		event.be_y[0] - yb,
// 		100.0
// 	);
// 	p_be = p_be.Unit() * be_momentum;

// 	// 4He momentum
// 	double he_momentum = MomentumFromKinetic(mass_4he, he_kinetic);
// 	// 4He momentum vector
// 	ROOT::Math::XYZVector p_he(
// 		event.he_x[1] - xb,
// 		event.he_y[1] - yb,
// 		100.0
// 	);
// 	p_he = p_he.Unit() * he_momentum;

// 	// excited 14C momentum vector
// 	ROOT::Math::XYZVector p_c = p_be + p_he;
// 	// excited 14C momentum
// 	double c_momentum = p_c.R();
// 	// excited 14C total energy
// 	double c_energy =
// 		(be_kinetic + mass_10be + excited_10be)
// 		+ (he_kinetic + mass_4he);
// 	// excited 14C mass
// 	double excited_c_mass = sqrt(
// 		pow(c_energy, 2.0) - pow(c_momentum, 2.0)
// 	);
// 	// excited energy of 14C
// 	double excited_14c = excited_c_mass - mass_14c;

// 	return excited_14c;
// }


class CostFunctor {
public:

	CostFunctor(
		const ThreeBodyInfoEvent &event,
		const int state
		// elc::DeltaEnergyCalculator *be_calculator,
		// elc::DeltaEnergyCalculator *he_calculator
	)
	: event_(event)
	, state_(state)
	// , be_calculator_(be_calculator)
	// , he_calculator_(he_calculator)
	{}


	bool operator()(
		const double * const csi_opt_param,
#ifdef POS_OPT
		const double * const pos_opt_param,
#endif
#ifdef PPAC_OPT
		const double * const ppac_param,
#endif
		double *residual
	) const {

#ifndef THRES_OPT
		double csi_param_full[3] = {csi_opt_param[0], csi_opt_param[1], 0.0};
#else
		double csi_param_full[3] = {
			csi_opt_param[0], csi_opt_param[1], csi_opt_param[2]
		};
#endif

#ifndef POS_OPT
		double pos_opt_param[2] = {0.0, 0.0};
#endif

// #ifndef PPAC_OPT
// 		double ppac_param[2] = {0.0, 0.0};
// #endif

		double d_kinetic, c_kinetic;

		double q = ThreeBodyProcess(
			event_,
			csi_param_full,
			pos_opt_param,
			d_kinetic, c_kinetic
		);

		double excited_energy = 0.0;
		if (state_ == 1) excited_energy = 3.368;
		else if (state_ == 2) excited_energy = 6.179;
		else if (state_ == 3) excited_energy = 6.179;
		// else if (state_ == 3) excited_energy = 7.542;

		residual[0] = q + (12.012512 + excited_energy);

		return true;
	}

private:
	ThreeBodyInfoEvent event_;
	int state_;
	// elc::DeltaEnergyCalculator *be_calculator_;
	// elc::DeltaEnergyCalculator *he_calculator_;
};


void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options]\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -c                Set calculate T0D2 flags.\n"
		"  -i                Iteration mode.\n";
}


/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] calculate flag of calculate T0D2
/// @param[out] iteration iteration mode
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	int &calculate,
	bool &iteration
) {
	// initialize
	help = false;
	calculate = 0;
	iteration = false;
	// start index of positional arugments
	int result = 0;
	for (result = 1; result < argc; ++result) {
		// assumed that all options have read
		if (argv[result][0] != '-') break;
		// short option contains only one letter
		if (argv[result][2] != 0) return -result;
		if (argv[result][1] == 'h') {
			help = true;
			return result;
		} else if (argv[result][1] == 'c') {
			// option of calculate T0D2 flag
			// get flag in next argument
			++result;
			// miss arguemnt behind option
			if (result == argc) return -argc;
			calculate = atoi(argv[result]);
		} else if (argv[result][1] == 'i') {
			// iteration mode
			iteration = true;
		}else {
			return -result;
		}
	}
	return result;
}


int main(int argc, char **argv) {
	google::InitGoogleLogging(argv[0]);

	// help flag
	bool help = false;
	// calculate flag
	int calculate_t0d2 = 0;
	// iteration mode
	bool iteration = false;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, calculate_t0d2, iteration);

	// need help
	if (help) {
		PrintUsage(argv[0]);
		return 0;
	}

	if (pos_start < 0) {
		if (-pos_start < argc) {
			std::cerr << "Error: Invaild option " << argv[-pos_start] << ".\n";
		} else {
			std::cerr << "Error: Option need parameter.\n";
		}
		PrintUsage(argv[0]);
		return -1;
	}

	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody.root", kGenerateDataPath, kInformationDir
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

	// resist strip file name
	TString resist_strip_file_name = TString::Format(
		"%s%sresist-strip.root", kGenerateDataPath, kShowDir
	);
	// add resist strip friend
	ipt->AddFriend("resist=tree", resist_strip_file_name);
	// most same state
	int most_same_state;
	// setup input branch
	ipt->SetBranchAddress("resist.same_state", &most_same_state);

	// output file
	TFile opf(
		TString::Format(
			"%s%sthreebody%s.root",
			kGenerateDataPath,
			kOptimizeDir,
			iteration ? "-1" : ""
		),
		"recreate"
	);
	// initial total Q value spectrum
	TH1F init_q_spectrum("hqi", "Q", 90, -23, -8);
	// initial Q value spectrum seperated by CsI index
	std::vector<TH1F> init_q_sep_spectrum;
	for (int i = 0; i < 12; ++i) {
		init_q_sep_spectrum.emplace_back(
			TString::Format("hqi%d", i), "Q",
			90, -23, -8
		);
	}
	// initial beam energy histogram
	TH1F init_hist_beam_energy("hbei", "beam energy", 100, 350, 400);
	// optimized total Q value spectrum
	TH1F opt_q_spectrum("hqo", "Q", 90, -23, -8);
	// optimized Q value spectrum seperated by CsI index
	std::vector<TH1F> opt_q_sep_spectrum;
	for (int i = 0; i < 12; ++i) {
		opt_q_sep_spectrum.emplace_back(
			TString::Format("hqo%d", i), "Q",
			90, -23, -8
		);
	}
	// optimized beam energy histogram
	TH1F opt_hist_beam_energy("hbeo", "beam energy", 100, 350, 400);
	// output tree
	TTree opt("tree", "threebody tree");
	// output data
	double be_kinetic;
	double he_kinetic;
	double d_kinetic;
	double c_kinetic;
	double q_value;
	int be_state;
	double c_excited;
	// setup output branches
	opt.Branch("csi_index", &event.csi_index, "ci/I");
	opt.Branch("be_kinetic", &be_kinetic, "bek/D");
	opt.Branch("he_kinetic", &he_kinetic, "hek/D");
	opt.Branch("d_kinetic", &d_kinetic, "dk/D");
	opt.Branch("c_kinetic", &c_kinetic, "ck/D");
	opt.Branch("q", &q_value, "q/D");
	opt.Branch("state", &be_state, "state/I");
	opt.Branch("excited", &c_excited, "ex/D");

	double csi_opt_param[12][3];
	double pos_opt_param[6][2];
	// double ppac_opt_param[2] = {0.0, 0.0};
	// double ppac_correct[2][3];

	// initialize
	for (int i = 0; i < 12; ++i) {
		for (int j = 0; j < 3; ++j) csi_opt_param[i][j] = csi_param[i][j];
	}
	for (int i = 0; i < 6; ++i) {
		pos_opt_param[i][0] = pos_opt_param[i][1] = 0.0;
	}
	// // read CsI calibration parameters from file
	// std::ifstream fin_csi(TString::Format(
	// 	"%s%scsi-calibrate%s.txt",
	// 	kGenerateDataPath,
	// 	kOptimizeDir,
	// 	iteration ? "-tb" : ""
	// ).Data());
	// for (int i = 0; i < 12; ++i) {
	// 	for (int j = 0; j < 3; ++j) fin_csi >> csi_param[i][j];
	// }
	// fin_csi.close();

	// // read TAF position parameters from file
	// std::ifstream fin_pos(TString::Format(
	// 	"%s%staf-position%s.txt",
	// 	kGenerateDataPath,
	// 	kOptimizeDir,
	// 	iteration ? "-tb" : ""
	// ).Data());
	// for (int i = 0; i < 6; ++i) {
	// 	for (int j = 0; j < 2; ++j) fin_pos >> pos_param[i][j];
	// }
	// fin_pos.close();

	// // read PPAC correct parameters from file
	// std::ifstream fin_ppac(TString::Format(
	// 	"%s%sppac-position%s.txt",
	// 	kGenerateDataPath,
	// 	kOptimizeDir,
	// 	iteration ? "-tb" : ""
	// ).Data());
	// for (int i = 0; i < 2; ++i) {
	// 	for (int j = 0; j < 3; ++j) fin_ppac >> ppac_correct[i][j];
	// }
	// fin_ppac.close();
	// ppac_param[0] = ppac_correct[0][0];
	// ppac_param[1] = ppac_correct[1][0];

	// energy calculators
	// elc::DeltaEnergyCalculator be_calculator("t0", "10Be");
	// elc::DeltaEnergyCalculator he_calculator("t0", "4He");


	int states_csi[12][5];
	for (int i = 0; i < 12; ++i) {
		for (int j = 0; j < 5; ++j) {
			states_csi[i][j] = 0;
		}
	}
	// std::vector<int> states;
	// loop to calculate Q value
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		// ignore events without PID
		if (event.taf_flag != 0) continue;
		// ignore events without PPAC tracking
		if (event.xppac_track[0] < 1 || event.xppac_track[1] < 1) continue;

		// q_value = ThreeBodyProcess(
		// 	event, (calculate_t0d2 & 0x1) != 0,
		// 	&be_calculator, &he_calculator,
		// 	csi_param[event.csi_index],
		// 	pos_param[event.csi_index/2],
		// 	ppac_param,
		// 	be_kinetic, he_kinetic, d_kinetic, c_kinetic
		// );

		init_q_spectrum.Fill(event.q);
		init_q_sep_spectrum[event.csi_index].Fill(event.q);
		init_hist_beam_energy.Fill(event.c14_kinetic);

		// get state
		// int state;
		// if (q_value > -13.5 && q_value < -10) state = 0;
		// else if (q_value > -16 && q_value < -14.5 ) state = 1;
		// else if (q_value > -19 && q_value < -16.5) state = 2;
		// else state = -1;

		// states.push_back(most_same_state);
		states_csi[event.csi_index][most_same_state+1]++;
	}

	std::cout << "  -1  0  1  2  3\n";
	for (int i = 0; i < 12; ++i) {
		std::cout << i << " ";
		for (int j = 0; j < 5; ++j) std::cout << states_csi[i][j] << " ";
		std::cout << "\n";
	}


	// begin numerical optimization
	ceres::Problem problem;
	// loop to add observations points to problem model
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		// ignore event without PID
		if (event.taf_flag != 0) continue;
		// ignore events without PPAC tracking
		if (event.xppac_track[0] < 1 || event.xppac_track[1] < 1) continue;

		// check state
		// if (states[entry] < 0) continue;
		if (most_same_state < 0) continue;

		ceres::CostFunction *cost_function;
		cost_function = new ceres::NumericDiffCostFunction<
			CostFunctor,
			ceres::CENTRAL,
			1
#ifdef THRES_OPT
			, 3
#else
			, 2
#endif
#ifdef POS_OPT
			, 2
#endif
#ifdef PPAC_OPT
			, 2
#endif
		> (
			new CostFunctor(
				event, most_same_state
			)
		);
		problem.AddResidualBlock(
			cost_function,
			nullptr,
			csi_opt_param[event.csi_index]
#ifdef POS_OPT
			, pos_opt_param[event.csi_index/2]
#endif
#ifdef PPAC_OPT
			, ppac_param
#endif
		);
	}

	// solve optimization problem
	ceres::Solver::Options options;
	options.max_num_iterations = 100;
	options.linear_solver_type = ceres::DENSE_QR;
	options.minimizer_progress_to_stdout = true;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	std::cout << summary.BriefReport() << "\n";

	// // get l0 from ppac parameters
	// ppac_correct[0][0] = ppac_param[0];
	// ppac_correct[1][0] = ppac_param[1];

	// PpacOffsetX(ppac_correct[0][0], ppac_correct[0][1], ppac_correct[0][2]);
	// PpacOffsetY(ppac_correct[1][0], ppac_correct[1][1], ppac_correct[1][2]);

	// print optimization results
	std::cout << "Optimization calibration results:\n";
	for (int i = 0; i < 12; ++i) {
		std::cout << i << " "
			<< csi_opt_param[i][0] << " "
			<< csi_opt_param[i][1] << " "
			<< csi_opt_param[i][2] << "\n";
	}
	std::cout << "Optimized position results:\n";
	for (int i = 0; i < 6; ++i) {
		std::cout << i << " "
			<< pos_opt_param[i][0] << " "
			<< pos_opt_param[i][1] << "\n";
	}
	std::cout << "Optimized PPAC position results:\n";
	for (int i = 0; i < 2; ++i) {
		std::cout << ppac_correct[i][0] << ", "
			<< ppac_correct[i][1] << ", "
			<< ppac_correct[i][2] << "\n";
	}

	// save calibration parameters
	std::ofstream fout_calibrate(TString::Format(
		"%s%scsi-calibrate-tb%s.txt",
		kGenerateDataPath,
		kOptimizeDir,
		iteration ? "1" : ""
	).Data());
	for (int i = 0; i < 12; ++i) {
		fout_calibrate
			<< csi_opt_param[i][0] << " "
			<< csi_opt_param[i][1] << " "
			<< csi_opt_param[i][2] << "\n";
	}
	fout_calibrate.close();

	// save position parameters
	std::ofstream fout_position(TString::Format(
		"%s%staf-position-tb%s.txt",
		kGenerateDataPath,
		kOptimizeDir,
		iteration ? "1" : ""
	).Data());
	for (int i = 0; i < 6; ++i) {
		fout_position << pos_opt_param[i][0] << " "
			<< pos_opt_param[i][1] << "\n";
	}
	fout_position.close();

	// // save PPAC position parameters
	// std::ofstream fout_ppac(TString::Format(
	// 	"%s%sppac-position-tb%s.txt",
	// 	kGenerateDataPath,
	// 	kOptimizeDir,
	// 	iteration ? "1" : ""
	// ).Data());
	// for (int i = 0; i < 2; ++i) {
	// 	fout_ppac << ppac_correct[i][0] << " "
	// 		<< ppac_correct[i][1] << " "
	// 		<< ppac_correct[i][2] << "\n";
	// }
	// fout_ppac.close();


	// loop to calculate Q value
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		// ignore event without PID
		if (event.taf_flag != 0) continue;
		// ignore events without PPAC tracking
		if (event.xppac_track[0] < 1 || event.xppac_track[1] < 1) continue;

		be_kinetic = event.t0_energy[0];
		he_kinetic = event.t0_energy[1];

		q_value = ThreeBodyProcess(
			event,
			csi_opt_param[event.csi_index],
			pos_opt_param[event.csi_index/2],
			d_kinetic, c_kinetic
		);

		// if (q_value > -13.5 && q_value < -10) be_state = 0;
		// else if (q_value > -16 && q_value < -14.5 ) be_state = 1;
		// else if (q_value > -19 && q_value < -16.5) be_state = 2;
		// else be_state = -1;

		// double be_excited = 0.0;
		// if (be_state == 1) be_excited = 3.368;
		// else if (be_state == 2) be_excited = 6.179;

		// c_excited = TwoBodyProcess(
		// 	event,
		// 	ppac_param,
		// 	be_excited
		// );

		opt_q_spectrum.Fill(q_value);
		opt_q_sep_spectrum[event.csi_index].Fill(q_value);
		opt_hist_beam_energy.Fill(c_kinetic);

		// fill to tree
		opt.Fill();
	}

	std::vector<int> calculate_flag;
	if (iteration) {
		TString last_file_name = TString::Format(
			"%s%sthreebody.root",
			kGenerateDataPath,
			kOptimizeDir
		);
		TFile last_file(last_file_name, "read");
		std::vector<int> *flag = (std::vector<int>*)last_file.Get("calculate");
		if (flag) {
			for (size_t i = 0; i < flag->size(); ++i) {
				calculate_flag.push_back(flag->at(i));
			}
		}
		last_file.Close();
	}
	calculate_flag.push_back(calculate_t0d2);

	// save histograms
	opf.cd();
	init_q_spectrum.Write();
	for (auto &hist : init_q_sep_spectrum) hist.Write();
	init_hist_beam_energy.Write();
	opt_q_spectrum.Write();
	for (auto &hist : opt_q_sep_spectrum) hist.Write();
	opt_hist_beam_energy.Write();
	// save tree
	opt.Write();
	// save calculate flag
	opf.WriteObject(&calculate_flag, "calculate");
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}