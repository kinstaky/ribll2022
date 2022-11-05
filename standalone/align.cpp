#include <iostream>

#include "include/alignment.h"

using namespace ribll;

/// @brief print usage of this program
///
/// @param[in] name program name
/// 
void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run\n"
		"  run                run number\n"
		"Options:\n"
		"  -h                 Print this help information.\n";
	return;
}



int main(int argc, char **argv) {
	if (argc != 2) {
		PrintUsage(argv[0]);
		return -1;
	}

	int run = atoi(argv[1]);

	Alignment align(run, 30, 2'000'000, -10'000'000'000, 10'000'000'000);
	align.SetVerbose();
	if (align.Align()) {
		std::cerr << "Error: Align failed.\n";
		return -1;
	}

	return 0;  
}