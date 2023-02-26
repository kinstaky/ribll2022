#include "include/detectors.h"

#include <iostream>

namespace ribll {

std::shared_ptr<Detector> CreateDetector(
	const std::string &name,
	unsigned int run
) {
	if (name == "tof") {
		return std::make_shared<Tof>(run);
	} else if (name == "t0d1") {
		return std::make_shared<T0d1>(run);
	} else if (name == "t0d2") {
		return std::make_shared<T0d2>(run);
	} else if (name == "t0d3") {
		return std::make_shared<T0d3>(run);
	} else {
		std::cerr << "Error: Create detector " << name << " failed.\n";
		return nullptr;
	}
}

}