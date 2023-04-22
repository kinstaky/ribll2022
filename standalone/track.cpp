#include <iostream>
#include <vector>
#include <string>

#include "include/telescopes.h"
#include "include/detector/ppac.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run telescope [tolerance]\n"
		"  run               Set run number.\n"
		"  telescope         Set telescope name.\n"
		"  tolerance         Set angle tolerance, default is 0.0.\n"
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

	if (pos_start + 1 >= argc) {
		// positional arguments less than 3
		std::cerr << "Error: Miss telescope argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}


	unsigned int run = atoi(argv[pos_start]);
	std::string name = argv[pos_start+1];
	double tolerance = 0.0;
	if (pos_start + 2 < argc) tolerance = atof(argv[pos_start+2]);

	if (name == "xppac" || name == "vppac") {
		Ppac ppac(run, name, tag);
		if (ppac.Track()) {
			std::cerr << "Error: Track " << name << " failed.\n";
			return -1;
		}
	} else {
		std::shared_ptr<Telescope> telescope =
			CreateTelescope(name, run, tag);
		if (!telescope) {
			std::cerr << "Error: Telescope " << name << " not found.\n";
			return -1;
		}
		if (telescope->Track(tolerance)) {
			std::cerr << "Error: Track " << name << " failed.\n";
			return -1;
		}
	}

	return 0;
}