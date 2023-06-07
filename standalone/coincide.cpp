#include <iostream>
#include <memory>

#include "include/channel.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run case\n"
		"  run               Set run number\n"
		"  case              Set the case:\n"
		"                      0 -- 8Be->4He+4He\n"
		"                      1 -- 12C->4He+4He+4He\n"
		"                      2 -- 14C->10Be+4He 2body\n"
		"                      3 -- 14C+2H->10Be+4He+2H 3body\n"
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
		channel = std::make_unique<C14ToBe10He4ThreeBodyChannel>(run);
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