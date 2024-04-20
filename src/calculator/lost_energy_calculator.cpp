#include "include/calculator/lost_energy_calculator.h"

#include <iostream>
#include <stdexcept>
#include <fstream>
#include <sstream>

#include "include/defs.h"

namespace ribll {

namespace elc {


LostEnergyCalculator::LostEnergyCalculator(
	const std::string &projectile
)
: projectile_(projectile) {

	input_file_ = nullptr;

	// input file
	TString input_file_name = TString::Format(
		"%s%s%s-Si-lost.root",
		kGenerateDataPath, kEnergyCalculateDir, projectile_.c_str()
	);
	input_file_ = new TFile(input_file_name, "read");
	if (!input_file_) {
		std::cerr << "Error: Open file "
			<< input_file_name << " failed.\n";
		throw std::runtime_error("input file not found.");
	}

	// read functions
	dedx_from_e_func_ = (TSpline3*)input_file_->Get("dedx_from_e");
}


LostEnergyCalculator::~LostEnergyCalculator() {
	if (input_file_) input_file_->Close();
}


int LostEnergyCalculator::Initialize(
	const std::vector<std::string> &projectiles
) {
	// data points
	int points = 0;
	// points' energy, x position
	double energy[264];
	// points' energy loss, y position
	double loss[264];
	for (const std::string &projectile : projectiles) {
		std::cout << "Energy lost " << projectile << "\n";

		// get projectile mass number
		std::stringstream ss;
		ss << projectile;
		int mass = 0;
		ss >> mass;

		// lise file to read
		TString lise_file_name = TString::Format(
			"%s%s%s-Si-lost.txt",
			kGenerateDataPath,
			kEnergyCalculateDir,
			projectile.c_str()
		);
		std::ifstream fin(lise_file_name.Data());
		if (!fin.good()) {
			std::cerr << "Error: Open file "
				<< lise_file_name << " failed.\n";
			return -1;
		}
		// read file
		std::string line;
		// read 3 useless lines
		for (int i = 0; i < 3; ++i) std::getline(fin, line);
		// read data
		points = 0;
		while (fin.good()) {
			std::getline(fin, line);
			std::stringstream ss;
			ss.str(line);
			double e, de;
			// get data from the second method
			ss >> e >> de >> e >> de;
			energy[points] = e * mass;
			loss[points] = de;
			++points;
			// stop at 500 MeV/u
			if (e >= 500) break;
		}
		fin.close();

		// output file name
		TString output_file_name = TString::Format(
			"%s%s%s-Si-lost.root",
			kGenerateDataPath, kEnergyCalculateDir, projectile.c_str()
		);
		// output file
		TFile opf(output_file_name, "recreate");

		double valbeg = (loss[1] - loss[0]) / (energy[1] - energy[0]);
		double valend = (loss[points-1] - loss[points-2])
			/ (energy[points-1] - energy[points-2]);
		// dE/dx from E function
		TSpline3 *func = new TSpline3(
			"dedx_from_e", energy, loss, points, "b1e1", valbeg, valend
		);
		// save
		func->Write("dedx_from_e");
		opf.Close();
	}
	return 0;
}


}	// elc

}	// ribll