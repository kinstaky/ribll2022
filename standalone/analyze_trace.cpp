#include <iostream>
#include <string>
#include <vector>

#include "include/detectors.h"
#include "include/telescopes.h"

using namespace ribll;


void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run detector\n"
		"  run               Set run number.\n"
		"  detector          Set detector name.\n"
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
	if (argc < 2) {
		PrintUsage(argv[0]);
		return -1;
	}

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

	if (pos_start+1 >= argc) {
		// positional arguments less than 2
		std::cerr << "Error: Miss detector argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}


	int run = atoi(argv[pos_start]);
	// list of detector names
	std::vector<std::string> detector_names;
	for (int i = pos_start+1; i < argc; ++i) {
		detector_names.push_back(std::string(argv[i]));
	}

	for (auto name : detector_names) {
		if (name.size() == 5 && name.substr(0, 4) == "tafd") {
			unsigned int index = name[4] - '0';
			if (index < 2 || index > 5) {
				std::cerr << "Error: Invalid detector " << name
					<< " for trace analysis.\n";
				continue;
			}
			Tafd tafd(run, index, "");
			if (tafd.AnalyzeTrace()) {
				std::cerr << "Error: Analyze trace of " << name
					<< " failed.\n";
				continue;
			}
		} else if (name.size() == 4 && name.substr(0, 3) == "taf") {
			unsigned int index = name[3] - '0';
			if (index > 5) {
				std::cerr << "Error: Invalid telescope " << name
					<< " for trace analysis.\n";
				continue;
			}
			Taf taf(run, index, "ta");
			if (taf.AnalyzeTrace()) {
				std::cerr << "Error: Analyze trace of " << name
					<< " failed.\n";
				continue;
			}
		}
	}
	return 0;
}