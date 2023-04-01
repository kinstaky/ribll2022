#include "include/telescope/telescope.h"

#include <iostream>

namespace ribll {

Telescope::Telescope(
	unsigned int run,
	const std::string &name,
	const std::string &tag
)
: run_(run)
, name_(name)
, tag_(tag) {

}


int Telescope::Track() {
	std::cerr << "Error: Not implemented yet.\n";
	return -1;
}

}		// namespace ribll