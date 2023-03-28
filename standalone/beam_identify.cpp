#include <iostream>

#include "include/detector/tof.h"

using namespace ribll;

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " run\n"
			<< "  run                   run number\n";
		return -1;
	}
	unsigned int run = atoi(argv[1]);

	Tof tof(run, "");
	if (tof.BeamIdentify()) {
		std::cerr << "Error: Identify beam falied.\n";
		return -1;
	}
	return 0;
}