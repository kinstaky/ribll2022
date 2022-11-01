#include <iostream>

#include "include/mapper.h"

using namespace ribll;


/// @brief print usage of this program
///
/// @param[in] name program name
/// 
void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run [crate]\n"
		"  run                run number\n"
		"  crate              Choose crate to mapping, 1, 2, 3 for xia and 4 for vme\n"
		"                       default(without this argument) is all.\n"
		"Options:\n"
		"  -h                 Print this help information.\n";
	return;
}


int main(int argc, char **argv) {
	if (argc < 2 || argc > 3) {
		PrintUsage(argv[0]);
		return -1;
	}

	int run = atoi(argv[1]);
	int crate = argc == 3 ? atoi(argv[2]) : 0;
	
	if (crate == 0 || crate == 1) {
		Crate1Mapper crate1(run);
		crate1.Mapping();
	}

	if (crate == 0 || crate == 2) {
		Crate2Mapper crate2(run);
		if (crate2.Mapping()) {
			std::cerr << "Error: mapping crate 2 failed.\n";
		}
	}

	// if (crate == 0 || crate == 3) {
	// 	Crate3Mapper crate3(run);
	// 	crate3.Map();
	// }

	if (crate == 0 || crate == 4) {
		Crate4Mapper crate4(run);
		crate4.Mapping();
	}
	return 0;
}