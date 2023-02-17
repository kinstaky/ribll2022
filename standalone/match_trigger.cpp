#include <iostream>
#include <string>

#include "include/detector/tof.h"

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

	if (detector_name == "tof") {
		Tof tof(run);
		if (tof.MatchTrigger(-1000, 1000)) {
			std::cerr << "Error: match trigger for "
				<< detector_name << " failed.\n";
			return -1;
		}
	} else if (detector_name == "t0d1") {
		std::cerr << "Error: t0d1 to do.\n";
		return -1;
	} else {
		std::cerr << "Error: match trigger in "
			<< detector_name << " failed.\n";
		return -1;
	}
	return 0;
}