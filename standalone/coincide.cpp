#include <iostream>
#include <memory>

#include "include/channel.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run case\n"
		"  run               Set run number\n"
		"  case              Set the case:\n"
		"                     -1 -- TAF particles\n"
		"                      0 -- 8Be->4He+4He\n"
		"                      1 -- 12C->4He+4He+4He\n"
		"                      2 -- 14C->10Be+4He 2body\n"
		"                      3 -- 14C+2H->10Be+4He+2H 3body\n"
		"                      4 -- 15C+1H->14C+2H\n"
		"                      5 -- 14C+1H->10Be+4He+1H 3body\n"
		"                      6 -- 14C->4He+4He+6He\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -r mass           Set recoil mass to choose from 1H, 2H or 3H.\n"
		"  -s                Use simulated data.\n";
}


/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] recoil_mass set recoil mass
/// @param[out] sim use simulated data
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	int &recoil_mass,
	bool &sim
) {
	// initialize
	help = false;
	recoil_mass = 0;
	sim = false;
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
		} else if (argv[result][1] == 'r') {
			// option of recoil mass
			// get mass in next argument
			++result;
			if (result == argc) return -argc;
			recoil_mass = atoi(argv[result]);
		} else if (argv[result][1] == 's') {
			sim = true;
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
	// recoil mass flag
	int recoil_mass = 0;
	// simulated data
	bool sim = false;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, recoil_mass, sim);
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
		std::cerr << "Error: Parameter case not found.\n";
		PrintUsage(argv[0]);
		return -1;
	}
	// run number
	unsigned int run = atoi(argv[pos_start]);
	int condition = atoi(argv[pos_start + 1]);

	// check run
	if (run == 628) {
		std::cerr << "Error: Bad run.\n";
		return -1;
	}

	if (condition == -1) {
		// merge TAF particles
		if (MergeTaf(run)) {
			std::cerr << "Error: Merge TAF particles failed.\n";
			return -1;
		} else {
			return 0;
		}
	}
	// coincide channels
	std::unique_ptr<Channel> channel;
	if (condition == 0) {
		// 8Be->4He+4He
		channel = std::make_unique<Be8ToTwoAlphaChannel>(run);
	} else if (condition == 1) {
		// 12C->4He+4He+4He
		channel = std::make_unique<C12ToThreeAlphaChannel>(run);
	} else if (condition == 2) {
		// 14C->10Be+4He 2body
		channel = std::make_unique<C14ToBe10He4TwoBodyChannel>(run);
	} else if (condition == 3) {
		// 14C+2H->10Be+4He+2H 3body
		// default
		if (recoil_mass == 0) recoil_mass = 2;
		channel = std::make_unique<C14ToBe10He4ThreeBodyChannel>(
			run, recoil_mass, sim
		);
	} else if (condition == 4) {
		// 15C+1H->14C+2H pd reaction
		channel = std::make_unique<C15pdChannel>(run);
	} else if (condition == 5) {
		// 14C+1H->10Be+4He+1H 3body
		channel = std::make_unique<C14ToBe10He4H1ThreeBodyChannel>(run);
	} else if (condition == 6) {
		// 14C->4He+4He+6He
		channel = std::make_unique<C14ToHe4He4He6Channel>(run);
	} else {
		std::cerr << "Error: Invalid case " << condition << ".\n";
		return -1;
	}
	if (channel->Coincide()) {
		std::cerr << "Error: Coincide in case " << condition << " failed.\n";
		return -1;
	}
	return 0;
}