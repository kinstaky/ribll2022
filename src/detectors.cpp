#include "include/detectors.h"

#include <iostream>

namespace ribll {

std::shared_ptr<Detector> CreateDetector(
	const std::string &name,
	unsigned int run
) {
	if (name == "tof") {
		return std::make_shared<Tof>(run);
	} else if (name == "vt") {
		return std::make_shared<VmeTrigger>(run);
	} else if (name == "xppac") {
		return std::make_shared<Ppac>(run, "xppac");
	} else if (name == "t0d1") {
		return std::make_shared<T0d1>(run);
	} else if (name == "t0d2") {
		return std::make_shared<T0d2>(run);
	} else if (name == "t0d3") {
		return std::make_shared<T0d3>(run);
	} else if (name.substr(0, 3) == "taf" && name.size() == 4) {
		unsigned int index = name[3] - '0';
		if (index <= 5) {
			return std::make_shared<Taf>(run, index);
		}
	} else if (name.substr(0, 3) == "tab" && name.size() == 4) {
		unsigned int index = name[3] - '0';
		if (index <= 5) {
			return std::make_shared<Tab>(run, index);
		}
	} else if (name == "tafcsi" || name == "tabcsi") {
		return std::make_shared<CircularCsi>(run, name);
	} else if (name == "t0csi" || name == "t1csi") {
		return std::make_shared<SquareCsi>(run, name);
	}

	std::cerr << "Error: Create detector " << name << " failed.\n";
	return nullptr;
}

}