#include "include/check/chain_check.h"

#include <iostream>
#include <string>

using namespace ribll;

int main(int argc, char **argv) {
	if (argc < 3) {
		std::cout << "Usage: " << argv[0] << " run type [addition]\n"
			<< "  run                run number\n"
			<< "  type               node type: map, align, match\n"
			<< "  addition           additional information for checking\n";
		return -1;
	}

	unsigned int run = atoi(argv[1]);
	std::string key = std::string(argv[2]);
	// add addition to key
	for (int i = 3; i < argc; ++i) {
		key += "-" + std::string(argv[i]);
	}

	if (!NodeCheck(key)) {
		std::cerr << "Error: Check type and addition invalid.\n";
		return -1;
	}

	if (ChainCheck(run, key)) {
		// print up to date
		std::cout << key << " is \033[0;32mup to date\033[0m.\n";
	}
	return 0;
}