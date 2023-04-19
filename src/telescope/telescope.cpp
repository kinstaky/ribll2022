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


int Telescope::Track(double) {
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


int Telescope::Particle() {
	std::cerr << "Error: Telescope::Particle not implemented yet.\n";
	return -1;
}


std::unique_ptr<TCutG> Telescope::ReadCut(
	const char *prefix,
	const char *particle
) const {
		// cut file name
	TString cut_file_name;
	cut_file_name.Form(
		"%s%scut/%s-%s-%s.txt",
		kGenerateDataPath,
		kParticleIdentifyDir,
		name_.c_str(),
		prefix,
		particle
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


int Telescope::ReadCalibrateParameters() {
	// parameters file name
	TString file_name;
	file_name.Form(
		"%s%s%s-calibration-param-%s%04u.txt",
		kGenerateDataPath,
		kCalibrationDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// parameters file
	std::ifstream fin(file_name.Data());
	// check file
	if (!fin.good()) {
		std::cerr << "Error: Read calibrate parameters from "
			<< file_name << " failed.\n";
		return -1;
	}
	// read parameters
	for (size_t i = 0; i < Layers(); ++i) {
		fin >> cali_params_[i*2] >> cali_params_[i*2+1];
	}
	// close file
	fin.close();
	return 0;
}


int Telescope::WriteCalibrateParameters() const {
	// parameters file name
	TString file_name;
	file_name.Form(
		"%s%s%s-calibration-param-%s%04u.txt",
		kGenerateDataPath,
		kCalibrationDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// parameters file
	std::ofstream fout(file_name.Data());
	// check file
	if (!fout.good()) {
		std::cerr << "Error: Open calibrate file "
			<< file_name << " failed.\n";
		return -1;
	}
	// write parameters
	for (size_t i = 0; i < Layers(); ++i) {
		fout << cali_params_[i*2] << " " << cali_params_[i*2+1] << "\n";
	}
	// close file
	fout.close();
	return 0;
}


}		// namespace ribll