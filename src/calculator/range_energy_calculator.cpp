#include "include/calculator/range_energy_calculator.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include "include/defs.h"

namespace ribll {

namespace elc {


RangeEnergyCalculator::RangeEnergyCalculator(
	const std::string &projectile,
	const std::string &material
)
: projectile_(projectile)
, material_(material) {

	// get projectile mass number
	if (projectile == "1H") mass_ = 1.0072764520;
	else if (projectile == "2H")  mass_ = 2.0135531980;
	else if (projectile == "3H")  mass_ = 3.0155007014;
	else if (projectile == "3He") mass_ = 3.0149321622;
	else if (projectile == "4He") mass_ = 4.0015060943;
	else if (projectile == "6He") mass_ = 6.0177887354;
	else if (projectile == "6Li") mass_ = 6.0134771477;
	else if (projectile == "7Li") mass_ = 7.0143576949;
	else if (projectile == "8Li") mass_ = 8.0208405004;
	else if (projectile == "7Be") mass_ = 7.0147343973;
	else if (projectile == "9Be") mass_ = 9.0099887420;
	else if (projectile == "10Be")mass_ = 10.0113403769;
	else if (projectile == "10B") mass_ = 10.0101939628;
	else if (projectile == "11B") mass_ = 11.0065622673;
	else if (projectile == "12C") mass_ = 11.9967085206;
	else if (projectile == "13C") mass_ = 13.0000633559;
	else if (projectile == "14C") mass_ = 13.9999505089;
	else if (projectile == "15C") mass_ = 15.0073077289;

	input_file_ = nullptr;
	range_energy_func_ = nullptr;
	energy_range_func_ = nullptr;

	// input file name
	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-%s.root",
		kGenerateDataPath,
		kEnergyCalculateDir,
		projectile_.c_str(),
		material_.c_str()
	);
	// input file
	input_file_ = new TFile(input_file_name, "read");
	range_energy_func_ = (TSpline3*)input_file_->Get("re");
	energy_range_func_ = (TSpline3*)input_file_->Get("er");
}


RangeEnergyCalculator::~RangeEnergyCalculator() {
	if (input_file_) input_file_->Close();
}


int RangeEnergyCalculator::Initialize(
	const std::vector<ProjectileMaterial> &list
) {
	// data
	// projectile mass number
	int mass = 0;
	// points in energy-range function
	int points = 0;
	// points of energy, in MeV
	double energy[256];
	// points of range, in um
	double range[256];

	for (const auto &[projectile, material] : list) {
		std::cout << "Range " << projectile << ", " << material << "\n";
		// get projectile mass number
		std::stringstream ss;
		ss << projectile;
		ss >> mass;

		// lise file to read
		TString lise_file_name;
		lise_file_name.Form(
			"%s%s%s-%s.txt",
			kGenerateDataPath,
			kEnergyCalculateDir,
			projectile.c_str(),
			material.c_str()
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
			double e, r;
			// get data from the second method
			ss >> e >> r >> e >> r;
			energy[points] = e * mass;
			range[points] = r;
			++points;
			// stop at 500 MeV/u
			if (e >= 500) break;
		}
		fin.close();

		// output root file name
		TString root_file_name;
		root_file_name.Form(
			"%s%s%s-%s.root",
			kGenerateDataPath,
			kEnergyCalculateDir,
			projectile.c_str(),
			material.c_str()
		);
		// output file
		TFile *opf = new TFile(root_file_name, "recreate");

		double valbeg = (energy[1] - energy[0]) / (range[1] - range[0]);
		double valend = energy[points-1] - energy[points-2];
		valend /= range[points-1] - range[points-2];
		// energy-range function, get energy from range
		TSpline3 *er_func = new TSpline3(
			"er", range, energy, points, "b1e1", valbeg, valend
		);

		valbeg = 1.0 / valbeg;
		valend = 1.0 / valend;
		// range-energy function, get range from energy
		TSpline3 *re_func = new TSpline3(
			"re", energy, range, points, "b1e1", valbeg, valend
		);

		// save functions and close file
		er_func->Write("er");
		re_func->Write("re");
		opf->Close();
	}
	return 0;
}

}	// namespace elc

}	// namespace ribll
