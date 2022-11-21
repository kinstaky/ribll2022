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

	for (const std::string &detector_name : detectors) {
		if (detector_name == "ppac") {
			PPAC ppac(run);
			if (ppac.Correlate()) {
				std::cerr << "Error: correlate ppac failed.\n";
			}
		} else {
			std::shared_ptr<Detector> detector = CreateDetector(run, detector_name);
			if (!detector) {
				std::cerr << "Warning: invalid detector " << detector_name << "\n";
				continue;
			}
			if (normalize) {
				if (detector->NormalCorrelate()) {
					std::cerr << "Error: normal correlate " << detector << " failed.\n";
					continue;
				}
			} else {
				if (detector->Correlate(135, 300)) {
					std::cerr << "Error: correlate " << detector << " failed.\n";
					continue;
				}
				if (merge && detector->Merge(0.02)) {
					std::cerr << "Error: merge " << detector << " failed.\n";
					continue;
				}
			}
		}
	}

    return 0;
}
