#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>

#include "include/event/ta_event.h"
#include "include/event/particle_type_event.h"
#include "include/calculator/range_energy_calculator.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run\n"
		"  run               Set run number.\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -t tag            Set trigger tag.\n";
}


/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] trigger_tag trigger tag get from arguments
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	std::string &trigger_tag
) {
	// initialize
	help = false;
	trigger_tag.clear();
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
		} else if (argv[result][1] == 't') {
			// option of trigger tag
			// get tag in next argument
			++result;
			// miss arguemnt behind option
			if (result == argc) return -argc;
			trigger_tag = argv[result];
		} else {
			return -result;
		}
	}
	return result;
}


int main(int argc, char **argv) {
	if (argc < 1) {
		PrintUsage(argv[0]);
		return -1;
	}

	// help flag
	bool help = false;
	// trigger tag
	std::string tag;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, tag);

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

	if (pos_start >= argc) {
		// positional arguments less than 2
		std::cerr << "Error: Miss TAFD-index argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}

	// run number
	int run = atoi(argv[pos_start]);

	// energy-range calculator
	elc::RangeEnergyCalculator calculator("2H", "Si");

	// output file name
	TString output_file_name = TString::Format(
		"%s%stafd-deduce-thick-%04d.root",
		kGenerateDataPath, kShowDir, run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// histograms
	std::vector<TH1F> hist_thicks;
	for (int i = 0; i < 12; ++i) {
		hist_thicks.emplace_back(
			TString::Format("ht%d", i), "thickness",
			100, 0, 500
		);
	}

	// loop TAFD index
	for (int index = 0; index < 6; ++index) {
		// input file name
		TString input_file_name = TString::Format(
			"%s%staf%d-telescope-%s%04d.root",
			kGenerateDataPath,
			kTelescopeDir,
			index,
			tag.empty() ? "" : (tag+"-").c_str(),
			run
		);
		// input file
		TFile ipf(input_file_name, "read");
		// input tree
		TTree *ipt = (TTree*)ipf.Get("tree");
		// input event
		TaEvent event;
		// setup input branches
		event.SetupInput(ipt);

		// input pid file name
		TString pid_file_name = TString::Format(
			"%s%staf%d-particle-type-%s%04d.root",
			kGenerateDataPath,
			kParticleIdentifyDir,
			index,
			tag.empty() ? "" : (tag+"-").c_str(),
			run
		);
		// add friend
		ipt->AddFriend("pid=tree", pid_file_name);
		// input pid event
		ParticleTypeEvent pid_event;
		// setup input branches
		pid_event.SetupInput(ipt, "pid.");

		// read group parameters
		double a0, a1, a2;
		std::vector<double> csi_a_param[3];
		std::vector<double> csi_b_param[3];
		// CSI-A parameter file name
		TString csi_a_file_name = TString::Format(
			"%s%scsi-group-parameters-%d.txt",
			kGenerateDataPath, kOptimizeDir, index*2
		);
		std::ifstream csi_a_file(csi_a_file_name.Data());
		if (!csi_a_file.good()) {
			std::cerr << "Error: Open CsI A file failed.\n";
			return -1;
		}
		while (csi_a_file.good()) {
			csi_a_file >> a0 >> a1 >> a2;
			csi_a_param[0].push_back(a0);
			csi_a_param[1].push_back(a1);
			csi_a_param[2].push_back(a2);
		}
		csi_a_file.close();
		// CSI-B parameter file name
		TString csi_b_file_name = TString::Format(
			"%s%scsi-group-parameters-%d.txt",
			kGenerateDataPath, kOptimizeDir, index*2+1
		);
		std::ifstream csi_b_file(csi_a_file_name.Data());
		if (!csi_b_file.good()) {
			std::cerr << "Error: Open CsI A file failed.\n";
			return -1;
		}
		while (csi_b_file.good()) {
			csi_b_file >> a0 >> a1 >> a2;
			csi_b_param[0].push_back(a0);
			csi_b_param[1].push_back(a1);
			csi_b_param[2].push_back(a2);
		}
		csi_b_file.close();

		int groups = csi_a_param[0].size();
		if (groups != int(csi_b_param[0].size())) {
			std::cerr << "Error: parameter groups of TAF"
				<< index << " is different.\n";
			return -1;
		}

		// total number of entries
		long long entries = ipt->GetEntries();
		// 1/100 of entries, for showing process
		long long entry100 = entries / 100 + 1;
		// show start
		printf("Deducing thickness of TAF%d   0%%", index);
		fflush(stdout);
		for (long long entry = 0; entry < entries; ++entry) {
			// show process
			if (entry % entry100 == 0) {
				printf("\b\b\b\b%3lld%%", entry / entry100);
				fflush(stdout);
			}
			// get data
			ipt->GetEntry(entry);

			// ignore empty events
			if (event.num != 1) continue;
			// only keep 2H event
			if (pid_event.mass[0] != 2) continue;

			// loop parameters
			for (int i = 0; i < groups; ++i) {
				if (event.flag[0] == 0x3) {
					double csi_energy = pow(
						(
							event.energy[0][1] - csi_a_param[2][i]
						) / csi_a_param[0][i],
						1.0 / csi_a_param[1][i]
					);
					double csi_range = calculator.Range(csi_energy);
					double full_range = calculator.Range(
						event.energy[0][0] + csi_energy
					);
					double thickness = full_range - csi_range;
					thickness *= cos(event.theta[0]);
					hist_thicks[index*2].Fill(thickness);
				} else if (event.flag[0] == 0x5) {
					double csi_energy = pow(
						(
							event.energy[0][1] - csi_b_param[2][i]
						) / csi_b_param[0][i],
						1.0 / csi_b_param[1][i]
					);
					double csi_range = calculator.Range(csi_energy);
					double full_range = calculator.Range(
						event.energy[0][0] + csi_energy
					);
					double thickness = full_range - csi_range;
					thickness *= cos(event.theta[0]);
					hist_thicks[index*2+1].Fill(thickness);
				}
			}
		}
		// show finish
		printf("\b\b\b\b100%%\n");

		// close files
		ipf.Close();
	}

	for (int i = 0; i < 12; ++i) {
		TF1 *fit_func = new TF1(
			TString::Format("f%d", i), "gaus", 100, 200
		);
		fit_func->SetParameter(1, 180);
		fit_func->SetParameter(2, 20);
		hist_thicks[i].Fit(fit_func, "RQ+");

		std::cout << i/2 << "AB"[i%2] << " "
			<< fit_func->GetParameter(1) << "\n";
	}

	opf.cd();
	// save histograms
	for (auto &hist : hist_thicks) hist.Write();
	// close output file
	opf.Close();

	return 0;
}