#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TString.h>

#include "include/event/t0_event.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run end_run\n"
		"  run               Set run number.\n"
		"  end_run           Set the last run to chain, included.\n"
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
	if (pos_start+1 >= argc) {
		// positional arguments less than 3
		std::cerr << "Error: Miss detector argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}
	// run number
	unsigned int run = atoi(argv[pos_start]);
	// run length
	unsigned int end_run = atoi(argv[pos_start+1]);

	// input t0 chain
	TChain chain("t0", "t0 telescope");
	for (unsigned int i = run; i <= end_run; ++i) {
		if (i == 628 || i == 630) continue;
		chain.AddFile(TString::Format(
			"%s%st0-telescope-%s%04u.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			tag.empty() ? "" : (tag+"-").c_str(),
			i
		));
	}
	// input t0 telescope event
	T0Event t0_event;
	// setup input branches
	t0_event.SetupInput(&chain);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%st0-pid-%s%04u-%04u.root",
		kGenerateDataPath,
		kShowDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run,
		end_run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// 2D histogram of d1d2 pid
	TH2F d1d2_pid(
		"pid12", "d1d2 #DeltaE-E pid",
		2500, 0, 50000, 3000, 0, 60000
	);
	// 2D histogram of d1d2 pid (stop in d2)
	TH2F d1d2_pid_stop(
		"pid12s", "d1d2 #DeltaE-E pid",
		2500, 0, 50000, 3000, 0, 60000
	);
	// 2D histogram of d2d3 pid
	TH2F d2d3_pid(
		"pid23", "d2d3 #DeltaE-E pid",
		1750, 0, 35000, 2000, 0, 40000
	);
	// 2D histogram of d2d3 pid stop in d3
	TH2F d2d3_pid_stop(
		"pid23s", "d2d3 #DeltaE-E pid",
		1750, 0, 35000, 2000, 0, 40000
	);
	// 2D histogram of d3s1
	TH2D d3s1_pid(
		"pid3s1", "d3s1 #DeltaE-E pid",
		3000, 0, 60000, 1750, 0, 35000
	);
	// 2D histogram of d3s1 stop in s1
	TH2D d3s1_pid_stop(
		"pid3s1s", "d3s1 #DeltaE-E pid",
		3000, 0, 60000, 1750, 0, 35000
	);
	// 2D histogram of s1s2
	TH2D s1s2_pid(
		"pids1s2", "s1s2 #DeltaE-E pid",
		2000, 0, 40000, 3000, 0, 60000
	);
	// 2D histogram of s1s2 stop in s2
	TH2D s1s2_pid_stop(
		"pids1s2s", "s1s2 #DeltaE-E pid",
		2000, 0, 40000, 3000, 0, 60000
	);
	// 2D histogram of s2s3
	TH2D s2s3_pid(
		"pids2s3", "s2s3 #DeltaE-E pid",
		2000, 0, 40000, 2000, 0, 40000
	);

	// total number of entries
	long long entries = chain.GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling T0 pid   0%%");
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
		if (t0_event.num != 1) continue;
		if (t0_event.flag[0] == 0x3) {
			d1d2_pid.Fill(
				t0_event.energy[0][1], t0_event.energy[0][0]
			);
			d1d2_pid_stop.Fill(
				t0_event.energy[0][1], t0_event.energy[0][0]
			);
		} else if (t0_event.flag[0] == 0x7) {
			d1d2_pid.Fill(t0_event.energy[0][1], t0_event.energy[0][0]);
			d2d3_pid.Fill(t0_event.energy[0][2], t0_event.energy[0][1]);
			if (t0_event.ssd_flag == 0x0) {
				d2d3_pid_stop.Fill(
					t0_event.energy[0][2], t0_event.energy[0][1]
				);
			} else if (t0_event.ssd_flag == 0x1) {
				d3s1_pid.Fill(
					t0_event.ssd_energy[0], t0_event.energy[0][2]
				);
				d3s1_pid_stop.Fill(
					t0_event.ssd_energy[0], t0_event.energy[0][2]
				);
			} else if (t0_event.ssd_flag == 0x3) {
				d3s1_pid.Fill(
					t0_event.ssd_energy[0], t0_event.energy[0][2]
				);
				s1s2_pid.Fill(
					t0_event.ssd_energy[1], t0_event.ssd_energy[0]
				);
				s1s2_pid_stop.Fill(
					t0_event.ssd_energy[1], t0_event.ssd_energy[0]
				);
			} else if (t0_event.ssd_flag == 0x7) {
				d3s1_pid.Fill(
					t0_event.ssd_energy[0], t0_event.energy[0][2]
				);
				s1s2_pid.Fill(
					t0_event.ssd_energy[1], t0_event.ssd_energy[0]
				);
				s2s3_pid.Fill(
					t0_event.ssd_energy[2], t0_event.ssd_energy[1]
				);
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save pid histograms
	d1d2_pid.Write();
	d1d2_pid_stop.Write();
	d2d3_pid.Write();
	d2d3_pid_stop.Write();
	d3s1_pid.Write();
	d3s1_pid_stop.Write();
	s1s2_pid.Write();
	s1s2_pid_stop.Write();
	s2s3_pid.Write();
	// close files
	opf.Close();
	return 0;
}