#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "include/event/dssd_event.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run detector\n"
		"  run               Set run number.\n"
		"  detector          Set detector name.\n"
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

	if (pos_start < 0) {
		if (-pos_start < argc) {
			std::cerr << "Error: Invaild option " << argv[-pos_start] << ".\n";
		} else {
			std::cerr << "Error: Option need parameter.\n";
		}
		PrintUsage(argv[0]);
		return -1;
	}

	if (pos_start+1 >= argc) {
		// positional arguments less than 2
		std::cerr << "Error: Miss detector argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}


	int run = atoi(argv[pos_start]);
	std::string detector_name(argv[pos_start+1]);

	// input file name
	TString input_file_name = TString::Format(
		"%s%s%s-merge-%s%04d.root",
		kGenerateDataPath,
		kMergeDir,
		detector_name.c_str(),
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	// input event
	DssdMergeEvent merge;
	// setup input branches
	merge.SetupInput(ipt);

	// normalize result friend
	// input normalize result file name
	TString normalize_result_file_name = TString::Format(
		"%s%s%s-result-%s%04d.root",
		kGenerateDataPath,
		kNormalizeDir,
		detector_name.c_str(),
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// add friend
	ipt->AddFriend("r=tree", normalize_result_file_name);
	// normalize result event
	DssdFundamentalEvent result;
	result.SetupInput(ipt, "r.");

	// output file name
	TString output_file_name = TString::Format(
		"%s%s%s-merge-counts-%s%04d.root",
		kGenerateDataPath,
		kShowDir,
		detector_name.c_str(),
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "merge counts");
	// output data
	int num;
	bool front_merge, back_merge;
	int front_strip, back_strip;
	double front_energy, back_energy;
	double front_ratio, back_ratio;
	long long merge_entry;
	// setup output branches
	opt.Branch("num", &num, "num/I");
	opt.Branch("front_merge", &front_merge, "fm/O");
	opt.Branch("back_merge", &back_merge, "bm/O");
	opt.Branch("front_strip", &front_strip, "fs/I");
	opt.Branch("back_strip", &back_strip, "bs/I");
	opt.Branch("front_energy", &front_energy, "fe/D");
	opt.Branch("back_energy", &back_energy, "be/D");
	opt.Branch("front_ratio", &front_ratio, "fr/D");
	opt.Branch("back_ratio", &back_ratio, "br/D");
	opt.Branch("entry", &merge_entry, "entry/L");

	// total entries of all events
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Counting   0%%");
	fflush(stdout);
	// loop events
	for (long long entry = 0; entry < entries; ++entry) {
		// show proces
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get data
		ipt->GetEntry(entry);

		num = merge.hit;
		for (int i = 0; i < merge.hit; ++i) {
			// count front side
			if (merge.merge_tag[i] == 0 || merge.merge_tag[i] == 1) {
				// front one strip case
				int index = -1;
				// search for front strip
				for (int j = 0; j < 8; ++j) {
					if ((merge.flag[i] & (1<<j)) != 0) {
						index = j;
						break;
					}
				}
				// not found, impossible
				if (index == -1) {
					std::cerr << "Error: Entry " << entry << "\n";
					return -1;
				}
				// fill front 1 strip information
				front_merge = false;
				front_strip = result.front_strip[index];
				front_energy = result.front_energy[index];
				front_ratio = 1.0;
			} else if (merge.merge_tag[i] == 2 || merge.merge_tag[i] == 3) {
				// front adjacent strip case
				int index0 = -1;
				int index1 = -1;
				// search for first strip
				for (int j = 0; j < 8; ++j) {
					if ((merge.flag[i] & (1<<j)) != 0) {
						index0 = j;
						break;
					}
				}
				// search for second strip
				for (int j = index0+1; j < 8; ++j) {
					if ((merge.flag[i] & (1<<j)) != 0) {
						index1 = j;
						break;
					}
				}
				// not found, impossible
				if (index0 == -1 || index1 == -1) {
					std::cerr << "Error: Entry " << entry << "\n";
					return -1;
				}
				// fill front adjacent strips information
				front_merge = true;
				front_strip =
					result.front_strip[index0] < result.front_strip[index1] ?
					result.front_strip[index0] : result.front_strip[index1];
				front_energy =
					result.front_energy[index0] + result.front_energy[index1];
				front_ratio = result.front_energy[index0] / front_energy;
				front_ratio = front_ratio < 0.5 ? 1 - front_ratio : front_ratio;
			}

			// count back side
			if (merge.merge_tag[i] == 0 || merge.merge_tag[i] == 2) {
				// back one strip case
				int index = -1;
				// search for back strip
				for (int j = 8; j < 16; ++j) {
					if ((merge.flag[i] & (1<<j)) != 0) {
						index = j;
						break;
					}
				}
				// not found, impossible
				if (index == -1) {
					std::cerr << "Error: Entry " << entry << "\n";
					return -1;
				}
				// correct
				index -= 8;
				// fill back 1 strip information
				back_merge = false;
				back_strip = result.back_strip[index];
				back_energy = result.back_energy[index];
				back_ratio = 1.0;
			} else if (merge.merge_tag[i] == 1 || merge.merge_tag[i] == 3) {
				// back adjacent strip case
				int index0 = -1;
				int index1 = -1;
				// search for first strip
				for (int j = 8; j < 16; ++j) {
					if ((merge.flag[i] & (1<<j)) != 0) {
						index0 = j;
						break;
					}
				}
				// search for second strip
				for (int j = index0+1; j < 16; ++j) {
					if ((merge.flag[i] & (1<<j)) != 0) {
						index1 = j;
						break;
					}
				}
				// not found, impossible
				if (index0 == -1 || index1 == -1) {
					std::cerr << "Error: Entry " << entry << "\n";
					return -1;
				}
				// correct
				index0 -= 8;
				index1 -= 8;
				// fill back adjacent strips information
				back_merge = true;
				back_strip =
					result.back_strip[index0] < result.back_strip[index1] ?
					result.back_strip[index0] : result.back_strip[index1];
				back_energy =
					result.back_energy[index0] + result.back_energy[index1];
				back_ratio = result.back_energy[index0] / back_energy;
				back_ratio = back_ratio < 0.5 ? 1 - back_ratio : back_ratio;
			}
			// entry
			merge_entry = entry;
			// fill data
			opt.Fill();
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();

	return 0;
}