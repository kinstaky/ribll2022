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
		"  -t tag            Set trigger tag.\n"
		"  -c                Use calibrated energy.\n"
		"  -o                Show hole area pid.\n"
		"  -d                Use double event.\n";
}

/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] trigger_tag trigger tag get from arguments
/// @param[out] calibration use calibration energy
/// @param[out] hole show hole area
/// @param[out] two use double events
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	std::string &trigger_tag,
	bool &calibration,
	bool &hole,
	bool &two
) {
	// initialize
	help = false;
	trigger_tag.clear();
	calibration = false;
	hole = false;
	two = false;
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
		} else if (argv[result][1] == 'c') {
			// option of calibration tag
			calibration = true;
		} else if (argv[result][1] == 'o') {
			// option of hole
			hole = true;
		} else if (argv[result][1] == 'd') {
			two = true;
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
	// calibration option
	bool calibration = false;
	// hole option
	bool hole = false;
	// double option
	bool two = false;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(
		argc, argv, help, tag, calibration, hole, two
	);
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
		if (i == 628) continue;
		chain.AddFile(TString::Format(
			"%s%st0-telescope-%s%04u.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			tag.empty() ? "" : (tag+"-").c_str(),
			i
		));
	}
	// TFile ipf(
	// 	TString::Format(
	// 		"%s%st0-telescope-%04u.root",
	// 		kGenerateDataPath,
	// 		kTelescopeDir,
	// 		run
	// 	),
	// 	"read"
	// );
	// TTree *ipt = (TTree*)ipf.Get("tree");
	// input t0 telescope event
	T0Event t0_event;
	// setup input branches
	t0_event.SetupInput(&chain);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%st0-pid-%s%s%04u-%04u.root",
		kGenerateDataPath,
		kShowDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		hole ? "hole-" : "",
		run,
		end_run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// 2D histogram of d1d2 pid
	TH2F d1d2_pid(
		"d1d2", "d1d2 #DeltaE-E pid",
		2500, 0, 50000, 3000, 0, 60000
	);
	// 2D histogram of d1d2 pid (stop in d2)
	TH2F d1d2_pid_stop(
		"d1d2s", "d1d2 #DeltaE-E pid",
		2500, 0, 50000, 3000, 0, 60000
	);
	// 2D histogram of d2d3 pid
	TH2F d2d3_pid(
		"d2d3", "d2d3 #DeltaE-E pid",
		1750, 0, 35000, 2000, 0, 40000
	);
	// 2D histogram of d2d3 pid stop in d3
	TH2F d2d3_pid_stop(
		"d2d3s", "d2d3 #DeltaE-E pid",
		1750, 0, 35000, 2000, 0, 40000
	);
	// 2D histogram of d3s1
	TH2D d3s1_pid(
		"d3s1", "d3s1 #DeltaE-E pid",
		3000, 0, 60000, 1750, 0, 35000
	);
	// 2D histogram of d3s1 stop in s1
	TH2D d3s1_pid_stop(
		"d3s1s", "d3s1 #DeltaE-E pid",
		3000, 0, 60000, 1750, 0, 35000
	);
	// 2D histogram of s1s2
	TH2D s1s2_pid(
		"s1s2", "s1s2 #DeltaE-E pid",
		2000, 0, 40000, 3000, 0, 60000
	);
	// 2D histogram of s1s2 stop in s2
	TH2D s1s2_pid_stop(
		"s1s2s", "s1s2 #DeltaE-E pid",
		2000, 0, 40000, 3000, 0, 60000
	);
	// 2D histogram of s2s3
	TH2D s2s3_pid(
		"s2s3", "s2s3 #DeltaE-E pid",
		2000, 0, 40000, 2000, 0, 40000
	);
	// 2D histogram of d1d3 pid
	TH2F d1d3_pid(
		"d1d3", "d1d3 #DeltaE-E pid",
		1750, 0, 35000, 750, 0, 15000
	);
	// 2D histogram of d1d3 pid (stop in d3)
	TH2F d1d3_pid_stop(
		"d1d3s", "d1d3 #DeltaE-E pid",
		1750, 0, 35000, 750, 0, 15000
	);

	// double cali_param[6][2];
	// for (size_t i = 0; i < 6; ++i) {
	// 	cali_param[i][0] = 0.0;
	// 	cali_param[i][1] = 1.0;
	// }
	// if (calibration) {
	// 	// parameters file
	// 	std::ifstream fin(TString::Format(
	// 		"%s%s%s-calibration-param%s-%04u.txt",
	// 		kGenerateDataPath,
	// 		kCalibrationDir,
	// 		name_.c_str(),
	// 		tag_.empty() ? "" : ("-"+tag_).c_str(),
	// 		run_
	// 	).Data());
	// 	// check file
	// 	if (!fin.good()) {
	// 		std::cerr << "Error: Read calibrate parameters from "
	// 			<< file_name << " failed.\n";
	// 		return -1;
	// 	}
	// 	// read parameters
	// 	for (size_t i = 0; i < Layers(); ++i) {
	// 		fin >> cali_params_[i*2] >> cali_params_[i*2+1];
	// 	}
	// 	// close file
	// 	fin.close();
	// 	return 0;
	// }

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

		if (two && t0_event.num != 2) continue;
		for (int i = 0; i < t0_event.num; ++i) {
			if (hole && !t0_event.hole[i]) continue;
			else if (t0_event.hole[i]) continue;
			if (t0_event.flag[i] == 0x3) {
				d1d2_pid.Fill(
					t0_event.energy[i][1], t0_event.energy[i][0]
				);
				d1d2_pid_stop.Fill(
					t0_event.energy[i][1], t0_event.energy[i][0]
				);
			} else if (t0_event.flag[i] == 0x7) {
				d1d2_pid.Fill(t0_event.energy[i][1], t0_event.energy[i][0]);
				d2d3_pid.Fill(t0_event.energy[i][2], t0_event.energy[i][1]);
				// d1d3_pid.Fill(t0_event.energy[i][2], t0_event.energy[i][0]);
				if (t0_event.ssd_flag == 0x0) {
					d2d3_pid_stop.Fill(
						t0_event.energy[i][2], t0_event.energy[i][1]
					);
					// d1d3_pid_stop.Fill(
					// 	t0_event.energy[i][2], t0_event.energy[i][0]
					// );
				} else if (t0_event.ssd_flag == 0x1) {
					d3s1_pid.Fill(
						t0_event.ssd_energy[0], t0_event.energy[i][2]
					);
					d3s1_pid_stop.Fill(
						t0_event.ssd_energy[0], t0_event.energy[i][2]
					);
				} else if (t0_event.ssd_flag == 0x3) {
					d3s1_pid.Fill(
						t0_event.ssd_energy[0], t0_event.energy[i][2]
					);
				} else if (t0_event.ssd_flag == 0x7) {
					d3s1_pid.Fill(
						t0_event.ssd_energy[0], t0_event.energy[0][2]
					);
				}
			} else if (t0_event.flag[i] == 0x5) {
				d1d3_pid.Fill(
					t0_event.energy[i][2], t0_event.energy[i][0]
				);
				if (t0_event.ssd_flag == 0x0) {
					d1d3_pid_stop.Fill(
						t0_event.energy[0][2], t0_event.energy[0][0]
					);
				}
			}
		}
		if (t0_event.ssd_flag == 0x3) {
			s1s2_pid.Fill(
				t0_event.ssd_energy[1], t0_event.ssd_energy[0]
			);
			s1s2_pid_stop.Fill(
				t0_event.ssd_energy[1], t0_event.ssd_energy[0]
			);
		} else if (t0_event.ssd_flag == 0x7) {
			s1s2_pid.Fill(
				t0_event.ssd_energy[1], t0_event.ssd_energy[0]
			);
			s2s3_pid.Fill(
				t0_event.ssd_energy[2], t0_event.ssd_energy[1]
			);
		}
		// if ((t0_event.flag[0] & 0x5) == 0x5) {
		// 	d1d3_pid.Fill(
		// 		t0_event.energy[0][2], t0_event.energy[0][0]
		// 	);
		// 	if (t0_event.ssd_flag == 0x0) {
		// 		d1d3_pid_stop.Fill(
		// 			t0_event.energy[0][2], t0_event.energy[0][0]
		// 		);
		// 	}
		// }
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
	d1d3_pid.Write();
	d1d3_pid_stop.Write();
	// close files
	opf.Close();
	return 0;
}