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

#include "include/event/generate_event.h"
#include "include/event/detect_event.h"

using namespace ribll;


class CostFunctor {
public:

	CostFunctor(
		const DetectEvent &detect,
		const double generate_point,
		const double be_distance,
		const double he_distance,
		const double d_distance,
		const double ppac_distance
	) {
		kinetic_[0] = detect.be_kinetic;
		kinetic_[1] = detect.he_kinetic;
		kinetic_[2] = detect.d_kinetic;
		kinetic_[3] = 385.0;
		distance_[0] = be_distance;
		distance_[1] = he_distance;
		distance_[2] = d_distance;
		distance_[3] = ppac_distance;
		reaction_point_ = generate_point;
	}


	bool operator()(const double * const param, double *residual) const {
		// calculate
		double numerator = 0.0;
		double denominator = 0.0;
		for (int i = 0; i < 4; ++i) {
			numerator += param[i] * sqrt(kinetic_[i]) * distance_[i];
			denominator += param[i] * sqrt(kinetic_[i]);
		}
		double reaction_point = numerator / denominator;
		residual[0] = reaction_point - reaction_point_;
		return true;
	}

private:
	// kinetic energy of each particle
	double kinetic_[4];
	// distance of each particle
	double distance_[4];
	// generated reaction point
	double reaction_point_;
};


void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options]\n"
		"Options:\n"
		"  -h                Print this help information.\n";
}


/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help
) {
	// initialize
	help = false;
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
		} else {
			return -result;
		}
	}
	return result;
}


int main(int argc, char **argv) {
	google::InitGoogleLogging(argv[0]);

	// help flag
	bool help = false;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help);
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

	// detected data file name
	TString detect_file_name = TString::Format(
		"%s%sdetect.root",
		kGenerateDataPath,
		kSimulateDir
	);
	// detected data file
	TFile ipf(detect_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< detect_file_name << " failed.\n";
		return -1;
	}
	// add generate friend
	ipt->AddFriend("gen=tree", TString::Format(
		"%s%sgenerate.root", kGenerateDataPath, kSimulateDir
	));
	// detect event
	DetectEvent detect;
	// generate event
	GenerateEvent generate;
	// setup input branches
	detect.SetupInput(ipt);
	generate.SetupInput(ipt, "gen.");

	// optimization problems
	ceres::Problem problems[6];
	// single ppac track parameters
	double spta_params[6][4];
	// initialize single ppac track parameters
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 2; ++j) {
			spta_params[2*i+j][0] = sqrt(10.0) / 100.0;
			spta_params[2*i+j][1] = sqrt(4.0) / 100.0;
			spta_params[2*i+j][2] = sqrt(2.0) / 135.0;
		}
		spta_params[2*i][3] = -sqrt(14.0) / ppac_xz[i];
		spta_params[2*i+1][3] = -sqrt(14.0) / ppac_yz[i];
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling optimizing data   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get data
		ipt->GetEntry(entry);
		if (detect.valid != 7) continue;
		for (int i = 0; i < 3; ++i) {
			ceres::CostFunction *cost_function_x
				= new ceres::NumericDiffCostFunction<
					CostFunctor, ceres::CENTRAL, 1, 4
				> (
					new CostFunctor(
						detect,
						generate.target_x,
						detect.t0x[0][0],
						detect.t0x[1][0],
						detect.tafx,
						detect.ppacx[i]
					)
				);
			problems[2*i].AddResidualBlock(
				cost_function_x, nullptr, spta_params[2*i]
			);
			ceres::CostFunction *cost_function_y
				= new ceres::NumericDiffCostFunction<
					CostFunctor, ceres::CENTRAL, 1, 4
				> (
					new CostFunctor(
						detect,
						generate.target_y,
						detect.t0y[0][0],
						detect.t0y[1][0],
						detect.tafy,
						detect.ppacy[i]
					)
				);
			problems[2*i+1].AddResidualBlock(
				cost_function_y, nullptr, spta_params[2*i+1]
			);
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// solve problems
	for (int i = 0; i < 6; ++i) {
		ceres::Solver::Options options;
		options.max_num_iterations = 100;
		options.linear_solver_type = ceres::DENSE_QR;
		options.minimizer_progress_to_stdout = true;
		ceres::Solver::Summary summary;
		ceres::Solve(options, problems+i, &summary);
		std::cout << summary.BriefReport() << "\n";
	}
	// print optimization results
	for (int i = 0; i < 6; ++i) {
		std::cout << "xy"[i%2] << i/3 << " "
			<< spta_params[i][0] << " "
			<< spta_params[i][1] << " "
			<< spta_params[i][2] << " "
			<< spta_params[i][3] << "\n";
	}

	// save optimization results
	std::ofstream fout(TString::Format(
		"%s%ssingle-ppac-track-a.txt",
		kGenerateDataPath, kOptimizeDir
	));
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 4; ++j) {
			fout << spta_params[i][j] << " \n"[j==3];
		}
	}
	fout.close();

	// output file name
	TString output_file_name = TString::Format(
		"%s%ssingle-ppac-track-ao.root",
		kGenerateDataPath, kSimulateDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// dX of generated and optimized data
	TH1F *hist_dx[3];
	for (int i = 0; i < 3; ++i) {
		hist_dx[i] = new TH1F(
			TString::Format("hdx%d", i),
			"dX of generarted and optimized data",
			100, -10, 10
		);
	}
	// dY of generated and optimized data
	TH1F *hist_dy[3];
	for (int i = 0; i < 3; ++i) {
		hist_dy[i] = new TH1F(
			TString::Format("hdy%d", i),
			"dY of generarted and optimized data",
			100, -10, 10
		);
	}
	// output tree
	TTree opt("tree", "optimized data");
	// output data
	double opt_tx[3], opt_ty[3];
	// setup output branches
	opt.Branch("ppac_tx", opt_tx, "sptx[3]/D");
	opt.Branch("ppac_ty", opt_ty, "spty[3]/D");
	opt.Branch("valid", &detect.valid, "vaild/I");

	// show start
	printf("Filling result   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get data
		ipt->GetEntry(entry);
		if (detect.valid != 7) {
			opt.Fill();
			continue;
		}
		for (int i = 0; i < 3; ++i) {
			double numerator_x = 0.0;
			numerator_x +=
				spta_params[2*i][0]
				* sqrt(detect.be_kinetic)
				* detect.t0x[0][0];
			numerator_x +=
				spta_params[2*i][1]
				* sqrt(detect.he_kinetic)
				* detect.t0x[1][0];
			numerator_x +=
				spta_params[2*i][2]
				* sqrt(detect.d_kinetic)
				* detect.tafx;
			numerator_x +=
				spta_params[2*i][3]
				* sqrt(385.0)
				* detect.ppacx[i];
			double denominator_x = 0.0;
			denominator_x += spta_params[2*i][0] * sqrt(detect.be_kinetic);
			denominator_x += spta_params[2*i][1] * sqrt(detect.he_kinetic);
			denominator_x += spta_params[2*i][2] * sqrt(detect.d_kinetic);
			denominator_x += spta_params[2*i][3] * sqrt(385.0);
			opt_tx[i] = numerator_x / denominator_x;

			double numerator_y = 0.0;
			numerator_y +=
				spta_params[2*i+1][0]
				* sqrt(detect.be_kinetic)
				* detect.t0y[0][0];
			numerator_y +=
				spta_params[2*i+1][1]
				* sqrt(detect.he_kinetic)
				* detect.t0y[1][0];
			numerator_y +=
				spta_params[2*i+1][2]
				* sqrt(detect.d_kinetic)
				* detect.tafy;
			numerator_y +=
				spta_params[2*i+1][3]
				* sqrt(385.0)
				* detect.ppacy[i];
			double denominator_y = 0.0;
			denominator_y += spta_params[2*i+1][0] * sqrt(detect.be_kinetic);
			denominator_y += spta_params[2*i+1][1] * sqrt(detect.he_kinetic);
			denominator_y += spta_params[2*i+1][2] * sqrt(detect.d_kinetic);
			denominator_y += spta_params[2*i+1][3] * sqrt(385.0);
			opt_ty[i] = numerator_y / denominator_y;
			// fill to histogram
			hist_dx[i]->Fill(opt_tx[i] - generate.target_x);
			hist_dy[i]->Fill(opt_ty[i] - generate.target_y);
		}
		// fill tree
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// write histograms
	for (int i = 0; i < 3; ++i) hist_dx[i]->Write();
	for (int i = 0; i < 3; ++i) hist_dy[i]->Write();
	// write tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}