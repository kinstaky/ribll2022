/*
 * This program maps single crate. To map all crates and align different
 * crates, see mapping.cpp
 */

#include <iostream>

#include "include/mapper.h"

using namespace ribll;


/// @brief print usage of this program
///
/// @param[in] name program name
/// 
void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run crate\n"
		"  run               run number\n"
		"  crate             Choose crate to mapping, 0, 1, 2 for xia and 3 for vme\n"
		"                      default(without this argument) is all.\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -n                Map without threshold.\n"
		"  -i                Map crate3 independent to XIA.\n";
	return;
}


/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] no_threshold map without threshold
/// @param[out] independet map crate3 independent to XIA
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	bool &no_threshold,
	bool &independent
) {
	// initialize
	help = false;
	no_threshold = false;
	independent = false;
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
		} else if (argv[result][1] == 'i') {
			independent = true;
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
	// threshold flag
	bool no_threshold = false;
	// independent flag
	bool independent = false;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(
		argc, argv,
		help, no_threshold, independent
	);

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

	if (pos_start + 1 >= argc) {
		// positional arguments less than 2
		std::cerr << "Error: Miss crate argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}

	int run = atoi(argv[pos_start]);
	int crate = atoi(argv[pos_start+1]);

	if (crate == 0) {
		Crate0Mapper crate0(run);
		if (crate0.Map(!no_threshold)) {
			std::cerr << "Error: Mapping crate 0 failed.\n";
		}
	} else if (crate == 1) {
		Crate1Mapper crate1(run);
		if (crate1.Map(!no_threshold)) {
			std::cerr << "Error: Mapping crate 1 failed.\n";
		}
	} else if (crate == 2) {
		Crate2Mapper crate2(run);
		if (crate2.Map(!no_threshold)) {
			std::cerr << "Error: Mapping crate 2 failed.\n";
		}
	} else if (crate == 3) {
		Crate3Mapper crate3(run);
		if (crate3.Map(independent)) {
			std::cerr << "Error: Mapping crate 3 failed.\n";
		}
	} else {
		std::cerr << "Error: Invalid crate number " << crate << "\n";
		return -1;
	}

	return 0;
}