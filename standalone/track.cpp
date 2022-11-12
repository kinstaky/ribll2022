#include <iostream>

#include "include/detector/ppac.h"

using namespace ribll;

int main(int argc, char **argv) {
	if (argc != 2) {
		std::cout << "Usage: " << argv[0] << " run\n";
		return 0;
	}
	int run = atoi(argv[1]);
	ribll::PPAC ppac(run);
	// if (ppac.Correlate()) {
	// 	std::cerr << "Error: correlate ppac failed.\n";
	// 	return -1;
	// }
	if (ppac.Tracking()) {
		std::cerr << "Error: tracking ppac failed.\n";
		return -1;
	}

	return 0;
}