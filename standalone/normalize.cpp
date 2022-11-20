#include <iostream>
#include <vector>
#include <string>

#include "include/detector/detector.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run length detector\n"
		"  run                         run number\n"
		"  length                      number of run to chain\n"
		"  detector                    detector name\n"
		"Options:\n"
		"  -h                          Print this help information.\n"
		"  -i                          Iteration mode.\n";		
}


int main(int argc, char **argv) {
	if (argc < 4) {
		PrintUsage(argv[0]);
		return -1;
	}

	bool iteration = false;

	int pos_arg_start = 1;
	// check parameters
	if (argv[1][0] == '-') {
		switch (argv[1][1]) {
			case 'h':
				PrintUsage(argv[0]);
				return 0;
			case 'i':
				iteration = true;
				++pos_arg_start;
				break;
			default:
				std::cerr << "Error: invalid option " << argv[1] << "\n";
		}
	}


	int run = atoi(argv[pos_arg_start]);
	int length = atoi(argv[pos_arg_start+1]);
	std::string detector_name = std::string(argv[pos_arg_start+2]);

	if (detector_name == "t0d1") {
		T0D1 t0d1(run, "t0d1", 135, 300);
		if (t0d1.Normalize(length, 29, 36, iteration)) {
			std::cerr << "Error: normalize " << detector_name << " failed.\n";
			return -1;
		}
	}
	return 0;
}