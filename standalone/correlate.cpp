#include <iostream>
#include <vector>
#include <string>

#include "include/detector/detector.h"
#include "include/detector/ppac.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run detector[...]\n"
		"  run                      run number.\n"
		"  detector                 detector to correlates: ppac\n"
		"\n"
		"Options:\n"
		"  -h                       Print this help information.\n"
		"  -n                       Normalize energy with parameters from txt file.\n"
		"  -m                       Merge adjacent events.\n";
}

int main(int argc, char **argv) {
    if (argc < 3) {
        std::cerr << "Error: parameters number error, 2 parameters is required.\n";
        PrintUsage(argv[0]);
        return -1;
    }

	int run = 0;
	bool normalize = false;
	bool merge = false;
	std::vector<std::string> detectors;
	int pos_arg_start = 1;
	
	// check
	if (argv[1][0] == '-') {
		switch (argv[1][1]) {
			case 'h':
				PrintUsage(argv[0]);
				return 0;
			case 'n':
				normalize = true;
				++pos_arg_start;
				break;
			case 'm':
				merge = true;
				++pos_arg_start;
				break;
			default:
				std::cerr << "Error: invalid option " << argv[1] << "\n";
				return -1;
		}
	}
	
	run = atoi(argv[pos_arg_start]);
	for (int i = pos_arg_start+1; i < argc; ++i) {
		detectors.push_back(argv[i]);
	}

	for (const std::string &detector : detectors) {
		if (detector == "ppac") {
			PPAC ppac(run);
			if (ppac.Correlate()) {
				std::cerr << "Error: correlate ppac failed.\n";
			}
		} else if (detector == "t0d1") {
			T0D1 t0d1(run, detector, 135, 300);
			if (normalize) {
				if (t0d1.NormalCorrelate()) {
					std::cerr << "Error: normal correlate " << detector << " failed.\n";
					return -1;
				}
			} else {
				if (t0d1.Correlate()) {
					std::cerr << "Error: correlate " << detector << " failed.\n";				
					return -1;
				}
				if (merge && t0d1.Merge()) {
					std::cerr << "Error: merge " << detector << " failed.\n";
					return -1;
				}
			}
		// } else if (detector == "t0d3") {
		// 	T0D3 t0d3(run, detector, 135, 300);
		// 	if (t0d3.Correlate()) {
		// 		std::cerr << "Error: correlate " << detector << " failed.\n";				
		// 	}
		// 	if (merge && t0d3.Merge()) {
		// 		std::cerr << "Error: merge " << detector << " failed.\n";
		// 	}
		}
	}

    return 0;
}
