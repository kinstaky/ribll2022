#include "include/calculator/delta_energy_calculator.h"

#include <iostream>
#include <map>
#include <stdexcept>

#include "include/defs.h"
#include "include/calculator/range_energy_calculator.h"

namespace ribll {

namespace elc {

std::map<std::string, int> layers{
	{"t0", 5}
};


DeltaEnergyCalculator::DeltaEnergyCalculator(
	const std::string &telescope,
	const std::string &projectile
)
: telescope_(telescope)
, projectile_(projectile) {

	input_file_ = nullptr;

	auto search = layers.find(telescope);
	if (search == layers.end()) {
		std::cerr << "Error: Telescope " << telescope
			<< " in delta calculation  not found.\n";
		throw std::runtime_error("telescope not found");
	}
	int layer = search->second;

	// input file name
	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-delta-%s.root",
		kGenerateDataPath,
		kEnergyCalculateDir,
		telescope_.c_str(),
		projectile_.c_str()
	);
	// input file
	input_file_ = new TFile(input_file_name, "read");
	if (!input_file_) {
		std::cerr << "Error: Open file " <<
			input_file_name << " failed.\n";
		throw std::runtime_error("input file not found.");
	}

	// read functions
	for (int i = 0; i < layer; ++i) {
		de_e_funcs_.push_back(
			(TSpline3*)input_file_->Get(TString::Format("de_e_%d", i))
		);
		e_de_funcs_.push_back(
			(TSpline3*)input_file_->Get(TString::Format("e_de_%d", i))
		);
	}
}


DeltaEnergyCalculator::~DeltaEnergyCalculator() {
	if (input_file_) input_file_->Close();
}


int DeltaEnergyCalculator::Initialize(
	const std::string &telescope,
	const std::vector<double> &thickness,
	const std::vector<std::string> &projectiles
) {

	for (const std::string &projectile : projectiles) {
		// output file name
		TString output_file_name;
		output_file_name.Form(
			"%s%s%s-delta-%s.root",
			kGenerateDataPath,
			kEnergyCalculateDir,
			telescope.c_str(),
			projectile.c_str()
		);
		// output file
		TFile output_file(output_file_name, "recreate");

		// range-energy calculator
		RangeEnergyCalculator calculator(projectile, "Si");

		// calculate step of total energy, in MeV
		double step = 0.1;
		// loop layers to calculate dE and E
		for (size_t i = 0; i < thickness.size()-1; ++i) {
			// minimum value of total energy
			double min_total_energy = calculator.Energy(thickness[i]);
			// maximum value of total energy
			double max_total_energy =
				calculator.Energy(thickness[i]+thickness[i+1]);
			// list of energy loss in the first layer
			std::vector<double> delta_energy_list;
			// list of energy loss in the second layer
			std::vector<double> energy_list;
			// loop the total energy and calculate energy loss in the
			// fist layer and the second layer
			for (
				double total_energy = min_total_energy;
				total_energy <= max_total_energy;
				total_energy += step
			) {
				// incident range under the total energy
				double range = calculator.Range(total_energy);
				// incident range in the second layer
				range -= thickness[i];
				// goto larger energy if can't pass the first layer
				if (range <= 0) continue;
				// pass the second layer, stop
				if (range > thickness[i+1]) break;
				// energy lost in the second layer
				double energy = calculator.Energy(range);
				if (energy < 0.0 || energy > total_energy) continue;
				energy_list.push_back(energy);
				// energy lost in the first layer
				double delta_energy = total_energy - energy;
				delta_energy_list.push_back(delta_energy);
			}

			size_t points = energy_list.size();
			double valbeg =
				delta_energy_list[1] - delta_energy_list[0];
			valbeg /= energy_list[1] - energy_list[0];
			double valend =
				delta_energy_list[points-1] - delta_energy_list[points-2];
			valend /= energy_list[points-1] - energy_list[points-2];

			// dE-E function, get dE from E
			TSpline3 *de_e_func = new TSpline3(
				TString::Format("de_e_%ld", i),
				energy_list.data(), delta_energy_list.data(), points,
				"b1e1",
				valbeg, valend
			);

			// reverse the list, TSpline3 requests increasing x position
			std::reverse(energy_list.begin(), energy_list.end());
			std::reverse(delta_energy_list.begin(), delta_energy_list.end());
			// calculate valbeg and valend
			double tmp = 1.0 / valbeg;
			valbeg = 1.0 / valend;
			valend = tmp;
			// E-dE funcion, get E from dE
			TSpline3 *e_de_func = new TSpline3(
				TString::Format("e_de_%ld", i),
				delta_energy_list.data(), energy_list.data(), points,
				"b1e1",
				valbeg, valend
			);

			// save functions
			output_file.cd();
			e_de_func->Write(TString::Format("e_de_%ld", i));
			de_e_func->Write(TString::Format("de_e_%ld", i));
		}
		output_file.Close();
	}
	return 0;
}

}	// namespace elc

}	// namespace ribll