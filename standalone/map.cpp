/*
 * This program maps all crates and align xia crate and vme crate.
 */

#include <iostream>

#include "include/alignment.h"
#include "include/mapper.h"

/// @brief print usage of this program
///
/// @param[in] name program name
///
void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run\n"
		"  run               run number\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -n                Map without threshold.\n";
	return;
}


/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] no_threshold map without threshold
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	bool &no_threshold
) {
	// initialize
	help = false;
	no_threshold = false;
	// start index of positional arugments
	int result = 0;
	for (result = 1; result < argc; ++result) {
		// assumed that all options have read
		if (argv[result][0] != '-') break;
		// short option contains only one letter
		if (argv[result][2] != 0) return -result;
		if (argv[result][1] == 'h') {
			help = true;
			return 0;
		} else if (argv[result][1] == 'n') {
			no_threshold = true;
		} else {
			return -result;
		}
	}
	return result;
}


int main(int argc, char **argv) {
	if (argc == 1) {
		PrintUsage(argv[0]);
		return -1;
	}

	// help flag
	bool help = false;
	// threshold flag
	bool no_threshold = false;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, no_threshold);

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

	int run = atoi(argv[pos_start]);

	ribll::Crate0Mapper crate0(run);
	if (crate0.Map(!no_threshold)) {
		std::cerr << "Error: Map crate 0 failed.\n";
		return -1;
	}

	ribll::Alignment align(run, 30, 2000000, -10000000000, 10000000000);
	align.SetVerbose(false);
	if (align.Align()) {
		std::cerr << "Error: Align failed.\n";
		return -1;
	}

	ribll::Crate1Mapper crate1(run);
	if (crate1.Map(!no_threshold)) {
		std::cerr << "Error: Map crate 1 failed.\n";
		return -1;
	}

	ribll::Crate2Mapper crate2(run);
	if (crate2.Map(!no_threshold)) {
		std::cerr << "Error: Map crate 2 failed.\n";
		return -1;
	}

	ribll::Crate3Mapper crate3(run);
	if (crate3.Map()) {
		std::cerr << "Error: Map crate 3 failed.\n";
		return -1;
	}

	return 0;
}