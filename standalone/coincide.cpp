#include <iostream>

#include "include/channel.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run particle..\n"
		"  run               Set run number\n"
		"  particle          Set the particle to catch,"
		"at least 2, at most 3.\n"
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
	if (argc < 4) {
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
	if (pos_start + 2 >= argc) {
		// positional arguments less than 3
		std::cerr << "Error: At least two particles.\n";
		PrintUsage(argv[0]);
		return -1;
	}
	if (pos_start + 4 < argc) {
		// positional arguments greater than 4
		std::cerr << "Error: At most three particles.\n";
		PrintUsage(argv[0]);
		return -1;
	}
	// run number
	int run = atoi(argv[pos_start]);
	std::vector<std::string> particles;
	for (int i = pos_start + 1; i < argc; ++i) {
		particles.push_back(argv[i]);
	}

	// check run
	if (run == 628) {
		std::cerr << "Error: Bad run.\n";
		return -1;
	}

	if (particles.size() == 2) {
		T0Channel channel_a(run, particles);
		if (channel_a.Coincide()) {
			std::cerr << "Error: Coincide T0 failed.\n";
			return -1;
		}
		T0TAFChannel channel_b(run, particles);
		if (channel_b.Coincide()) {
			std::cerr << "Error: Coincide T0 TAF failed.\n";
			return -1;
		}
	} else if (particles.size() == 3) {
		T0TAFChannel channel(run, particles);
		if (channel.Coincide()) {
			std::cerr << "Error: Coincide T0 TAF failed.\n";
			return -1;
		}
	} else {
		std::cerr << "Error: Should not be here.\n";
		return -1;
	}
	return 0;
}