#include "include/telescope/telescope.h"

#include <iostream>

#include <TMath.h>

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
	std::cerr << "Error: Telescope::Track not implemented yet.\n";
	return -1;
}


int Telescope::ParticleIdentify() {
	std::cerr << "Error: Telescope::ParticleIdentify not implemented yet.\n";
	return -1;
}

int Telescope::Calibrate() {
	std::cerr << "Error: Telescope::Calibrate not implemented yet.\n";
	return -1;
}


int Telescope::AlphaCalibrate() {
	std::cerr << "Error: Telescope::AlphaCalibrate not implemented yet.\n";
	return -1;
}


int Telescope::CsiCalibrate() {
	std::cerr << "Error: Telescope::CsiCalibrate not implemented yet.\n";
	return -1;
}

}		// namespace ribll