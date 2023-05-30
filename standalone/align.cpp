#include <iostream>
#include <map>

#include "include/alignment.h"

using namespace ribll;

/// @brief print usage of this program
///
/// @param[in] name program name
///
void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run\n"
		"  run                Set run number.\n"
		"Options:\n"
		"  -h                 Print this help information.\n";
	return;
}

const std::map<int, int> group_map = {
	{624, 90}
};

int main(int argc, char **argv) {
	if (argc != 2) {
		PrintUsage(argv[0]);
		return -1;
	}

	// run number
	int run = atoi(argv[1]);

	int group = 100;
	auto search = group_map.find(run);
	if (search != group_map.end()) {
		group = search->second;
	}

	Alignment align(
		run, group, 10'000'000,
		-10'000'000'000, 10'000'000'000
	);
	align.SetVerbose(false);
	if (align.Align()) {
		std::cerr << "Error: Align failed.\n";
		return -1;
	}

	return 0;
}