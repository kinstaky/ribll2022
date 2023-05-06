#include <iostream>
#include <string>
#include <vector>

#include "include/mapper/trace_mapper.h"

using namespace ribll;


void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run crate\n"
		"  run               Set run number.\n"
		"  crate             Set crate id.\n"
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


	unsigned int run = atoi(argv[pos_start]);
	unsigned int crate = atoi(argv[pos_start+1]);

	TraceMapper mapper(run, crate);
	if (mapper.Map()) {
		std::cerr << "Error: Map trace of crate " << crate << " failed.\n";
		return -1;
	}

	return 0;
}