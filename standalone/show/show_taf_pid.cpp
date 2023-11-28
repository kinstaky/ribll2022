#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TString.h>

#include "include/event/ta_event.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run end_run index\n"
		"  run               Set run number.\n"
		"  end_run           Set the last run to chain, included.\n"
		"  index             Set the index of the TA telescope, 0-5.\n"
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
	if (argc < 3) {
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
	// invalid arguments
	if (pos_start < 0) {
		if (-pos_start < argc) {
			std::cerr << "Error: Invaild option " << argv[-pos_start] << ".\n";
		} else {
			std::cerr << "Error: Option need parameter.\n";
		}
		PrintUsage(argv[0]);
		return -1;
	}
	// check number of positional arguments
	if (pos_start+2 >= argc) {
		// positional arguments less than 3
		std::cerr << "Error: Miss detector argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}
	// run number
	unsigned int run = atoi(argv[pos_start]);
	// run length
	unsigned int end_run = atoi(argv[pos_start+1]);
	// index
	int index = atoi(argv[pos_start+2]);

	// input t0 chain
	TChain chain("ta", "ta telescope");
	for (unsigned int i = run; i <= end_run; ++i) {
		if (i == 628) continue;
		chain.AddFile(TString::Format(
			"%s%staf%d-telescope-%s%04u.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			index,
			tag.empty() ? "" : (tag+"-").c_str(),
			i
		));
	}
	// input t0 telescope event
	TaEvent ta_event;
	// setup input branches
	ta_event.SetupInput(&chain);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%staf%d-pid-%s%04u-%04u.root",
		kGenerateDataPath,
		kShowDir,
		index,
		tag.empty() ? "" : (tag+"-").c_str(),
		run,
		end_run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// 2D histogram of total pid
	TH2F total_pid(
		"pid", "TAFD-CsI(Tl) #DeltaE-E pid",
		3000, 0, 30000, 2000, 0, 20.0
	);
	// 2D histogram of pid under theta 0.8
	TH2F pid_l8(
		"pidl", "TAFD-CsI(Tl) #DeltaE-E pid #theta<0.8",
		3000, 0, 30000, 2000, 0, 20.0
	);
	// 2D histogram of pid above theta 0.8
	TH2F pid_g8(
		"pidg", "TAFD-CsI(Tl) #DeltaE-E pid #theta>0.8",
		3000, 0, 30000, 2000, 0, 20.0
	);
	// 2D histogram of pid under theta 0.8 and CsI-A
	TH2F pid_l8_a(
		"pidla", "TAFD-CsI(Tl)-A #DeltaE-E pid #theta<0.8",
		3000, 0, 30000, 2000, 0, 20.0
	);
	// 2D histogram of pid above theta 0.8 and CsI-A
	TH2F pid_g8_a(
		"pidga", "TAFD-CsI(Tl)-A #DeltaE-E pid #theta>0.8",
		3000, 0, 30000, 2000, 0, 20.0
	);
	// 2D histogram of pid under theta 0.8 and CsI-B
	TH2F pid_l8_b(
		"pidlb", "TAFD-CsI(Tl)-B #DeltaE-E pid #theta<0.8",
		3000, 0, 30000, 2000, 0, 20.0
	);
	// 2D histogram of pid above theta 0.8 and CsI-B
	TH2F pid_g8_b(
		"pidgb", "TAFD-CsI(Tl)-B #DeltaE-E pid #theta>0.8",
		3000, 0, 30000, 2000, 0, 20.0
	);
	// 2D dE-E PID histogram of CsI A with thick correctness
	TH2F pid_a_thick(
		"pidat", "TAFD-CsI(Tl)-A #DeltaE-E pid thick correctness",
		3000, 0, 30000, 2000, 0, 20.0
	);
	// 2D dE-E PID histogram of CsI B with thick correctness
	TH2F pid_b_thick(
		"pidbt", "TAFD-CsI(Tl)-B #DeltaE-E pid thick correctness",
		3000, 0, 30000, 2000, 0, 20.0
	);

	// 2D calibrated dE-E PID histogram of CsI A
	TH2F pid_a_calibrated(
		"pidac", "TAFD-CsI(Tl)-A #DeltaE-E calibrated",
		3000, 0, 100, 2000, 0, 20.0
	);
	// 2D calibrated dE-E PID histogram of CsI A
	TH2F pid_g8_a_calibrated(
		"pidgac", "TAFD-CsI(Tl)-A #DeltaE-E calibrated",
		3000, 0, 100, 2000, 0, 20.0
	);

	// total number of entries
	long long entries = chain.GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling TAF%d pid   0%%", index);
	fflush(stdout);
	// loop to fill pid histogram
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		chain.GetEntry(entry);
		if (ta_event.num != 1) continue;
		double &de = ta_event.energy[0][0];
		double &e = ta_event.energy[0][1];
		// delta energy after thick correct
		double cde = de * cos(ta_event.theta[0]);
		if (ta_event.flag[0] == 0x3) {
			total_pid.Fill(e, de);
			if (ta_event.theta[0] < 0.8) {
				pid_l8.Fill(e, de);
				pid_l8_a.Fill(e, de);
			} else {
				pid_g8.Fill(e, de);
				pid_g8_a.Fill(e, de);
				pid_g8_a_calibrated.Fill(pow((e+230.988)/245.278, 1.0/0.945929), de);
			}
			pid_a_thick.Fill(e+(de-cde)*150.0, cde);
			pid_a_calibrated.Fill(pow((e+230.988)/245.278, 1.0/0.945929), de);
		} else if (ta_event.flag[0] == 0x5) {
			if (ta_event.theta[0] < 0.8) {
				pid_l8.Fill(e, de);
				pid_l8_b.Fill(e, de);
			} else {
				pid_g8.Fill(e, de);
				pid_g8_b.Fill(e, de);
			}
			pid_b_thick.Fill(e+(de-cde)*150.0, cde);
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save pid histograms
	total_pid.Write();
	pid_l8.Write();
	pid_g8.Write();
	pid_l8_a.Write();
	pid_g8_a.Write();
	pid_l8_b.Write();
	pid_g8_b.Write();
	pid_a_thick.Write();
	pid_b_thick.Write();
	pid_a_calibrated.Write();
	pid_g8_a_calibrated.Write();
	// close files
	opf.Close();
	return 0;
}