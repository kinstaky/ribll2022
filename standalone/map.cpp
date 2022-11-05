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

	ribll::Crate1Mapper crate1(run);
	if (crate1.Map()) {
		std::cerr << "Error: Map crate 1 failed.\n";
		return -1;
	}

	ribll::Alignment align(run, 30, 2'000'000, -10'000'000'000, 10'000'000'000);
	align.SetVerbose(false);
	if (align.Align()) {
		std::cerr << "Error: Align failed.\n";
		return -1;
	}

	ribll::Crate2Mapper crate2(run);
	if (crate2.Map()) {
		std::cerr << "Error: Map crate 2 failed.\n";
		return -1;
	}

	ribll::Crate3Mapper crate3(run);
	if (crate3.Map()) {
		std::cerr << "Error: Map crate 3 failed.\n";
		return -1;
	}

	ribll::Crate4Mapper crate4(run);
	if (crate4.Map()) {
		std::cerr << "Error: Map crate 4 failed.\n";
		return -1;
	}

	return 0;
}