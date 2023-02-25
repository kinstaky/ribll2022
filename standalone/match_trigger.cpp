#include <iostream>
#include <string>

#include "include/detectors.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " run detector\n"
		<< "  run                       run number\n"
		<< "  detector                  detector name\n";
}

int main(int argc, char **argv) {
	// check parameters
	if (argc != 3) {
		PrintUsage(argv[0]);
		return -1;
	}

	int run = atoi(argv[1]);
	std::string detector_name = std::string(argv[2]);

	std::shared_ptr<Detector> detector = CreateDetector(detector_name, run);
	if (!detector) return -1;

	if (detector->MatchTrigger(-1000, 1000)) {
		std::cerr << "Error: match trigger for "
			<< detector_name << " failed.\n";
		return -1;
	}
	return 0;
}