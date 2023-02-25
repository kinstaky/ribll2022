#include "include/detectors.h"

#include <iostream>

namespace ribll {

std::shared_ptr<Detector> CreateDetector(
	const std::string &name,
	unsigned int run
) {
	if (name == "tof") {
		return std::make_shared<Tof>(run);
	} else {
		std::cerr << "Error: Create detector " << name << " failed.\n";
		return nullptr;
	}
}

}