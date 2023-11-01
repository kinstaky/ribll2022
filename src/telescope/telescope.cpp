#include "include/telescope/telescope.h"

#include <iostream>
#include <fstream>

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

int Telescope::Calibrate(unsigned int) {
	std::cerr << "Error: Telescope::Calibrate not implemented yet.\n";
	return -1;
}


int Telescope::AlphaCalibrate() {
	std::cerr << "Error: Telescope::AlphaCalibrate not implemented yet.\n";
	return -1;
}


int Telescope::CsiCalibrate(unsigned int) {
	std::cerr << "Error: Telescope::CsiCalibrate not implemented yet.\n";
	return -1;
}


int Telescope::Rebuild() {
	std::cerr << "Error: Telescope::Particle not implemented yet.\n";
	return -1;
}


std::unique_ptr<TCutG> Telescope::ReadCut(const char *name) const {
	// cut file name
	TString cut_file_name;
	cut_file_name.Form(
		"%s%scut/%s.txt",
		kGenerateDataPath,
		kParticleIdentifyDir,
		name
	);
	// open cut file to read points
	std::ifstream fin(cut_file_name.Data());
	if (!fin.good()) {
		std::cerr << "Error: Open file "
			<< cut_file_name << " failed.\n";
		return nullptr;
	}
	// result
	std::unique_ptr<TCutG> result = std::make_unique<TCutG>();
	// point index
	int point;
	// point positions
	double x, y;
	// loop to read points
	while (fin.good()) {
		fin >> point >> x >> y;
		result->SetPoint(point, x, y);
	}
	// close file
	fin.close();
	return result;
}


int Telescope::CalibrateResult() {
	std::cerr << "Error: Telescope::CalibrateResult is not implemented yet.\n";
	return -1;
}


}		// namespace ribll