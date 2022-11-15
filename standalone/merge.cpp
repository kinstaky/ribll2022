#include <iostream>
#include <string>

#include "include/detector/detector.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " run detector\n"
		"  run                          run number\n"
		"  detector                     detector name\n";
}

int main(int argc, char **argv) {
	if (argc != 3) {
		PrintUsage(argv[0]);
		return -1;
	}

	int run = atoi(argv[1]);
	std::string detector_name = std::string(argv[2]);

	if (detector_name == "t0d1") {
		T0D1 t0d1(run, "t0d1", 135, 300);
		if (t0d1.Merge()) {
			std::cerr << "Error: merge adjacent strip events in "
				<< detector_name << " failed.\n";
			return -1;
		}
	}
	return 0;
}