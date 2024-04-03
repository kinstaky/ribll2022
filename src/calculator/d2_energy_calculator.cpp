#include "include/calculator/d2_energy_calculator.h"

#include <iostream>

#include "include/defs.h"
#include "include/calculator/range_energy_calculator.h"

namespace ribll {

namespace elc {

D2EnergyCalculator::D2EnergyCalculator(const std::string &projectile)
: si_calculator_(projectile, "Si") {
	// input file name
	TString input_file_name;
	input_file_name.Form(
		"%s%st0-d2-%s.root",
		kGenerateDataPath,
		kEnergyCalculateDir,
		projectile.c_str()
	);
	// input file
	input_file_ = new TFile(input_file_name, "read");
	if (!input_file_) {
		std::cerr << "Error: Open file " <<
			input_file_name << " failed.\n";
		throw std::runtime_error("input file not found.");
	}

	// read functions
	for (int i = 0; i < 2; ++i) {
		de_e_stop_funcs_.push_back(
			(TSpline3*)input_file_->Get(TString::Format("de_e_stop_%d", i))
		);
		de_e_pass_funcs_.push_back(
			(TSpline3*)input_file_->Get(TString::Format("de_e_pass_%d", i))
		);
		e_de_stop_funcs_.push_back(
			(TSpline3*)input_file_->Get(TString::Format("e_de_stop_%d", i))
		);
		e_de_pass_funcs_.push_back(
			(TSpline3*)input_file_->Get(TString::Format("e_de_pass_%d", i))
		);
	}
}


D2EnergyCalculator::~D2EnergyCalculator() {
	if (input_file_) input_file_->Close();
}


int D2EnergyCalculator::Initialize(
	const std::vector<double> &thickness,
	const std::vector<std::string> &projectiles
) {
	for (const std::string &projectile : projectiles) {
		std::cout << "D2 " << projectile << "\n";
		// output file name
		TString output_file_name = TString::Format(
			"%s%st0-d2-%s.root",
			kGenerateDataPath,
			kEnergyCalculateDir,
			projectile.c_str()
		);
		// output file
		TFile output_file(output_file_name, "recreate");

		// range-energy calculator
		RangeEnergyCalculator calculator(projectile, "Si");

		std::vector<double> d1d2_stop_de_list;
		std::vector<double> d1d2_stop_e_list;
		std::vector<double> d1d2_pass_de_list;
		std::vector<double> d1d2_pass_e_list;
		std::vector<double> d2d3_stop_de_list;
		std::vector<double> d2d3_stop_e_list;
		std::vector<double> d2d3_pass_de_list;
		std::vector<double> d2d3_pass_e_list;

		// calculate step of total energy, in MeV
		double step = 0.1;
		double max_energy = 400.0;
		double min_energy = calculator.Energy(thickness[0]);
		for (
			double full_energy = min_energy;
			full_energy < max_energy;
			full_energy += step
		) {
			double range = calculator.Range(full_energy);
			if (
				range > thickness[0]
				&& range <= thickness[0]+thickness[1]
			) {
				double d2_energy = calculator.Energy(range-thickness[0]);
				double d1_energy = full_energy - d2_energy;
				d1d2_stop_de_list.push_back(d1_energy);
				d1d2_stop_e_list.push_back(d2_energy);
			} else if (
				range > thickness[0]+thickness[1]
				&& range <= thickness[0]+thickness[1]+thickness[2]
			) {
				double d3_energy =
					calculator.Energy(range-thickness[0]-thickness[1]);
				double d2d3_energy = calculator.Energy(range-thickness[0]);
				double d2_energy = d2d3_energy - d3_energy;
				double d1_energy = full_energy - d2d3_energy;
				d1d2_pass_de_list.push_back(d1_energy);
				d1d2_pass_e_list.push_back(d2_energy);
				d2d3_stop_de_list.push_back(d2_energy);
				d2d3_stop_e_list.push_back(d3_energy);
			} else if (
				range > thickness[0]+thickness[1]+thickness[2]
			) {
				double ssd_energy = calculator.Energy(
					range-thickness[0]-thickness[1]-thickness[2]
				);
				double d3ssd_energy = calculator.Energy(
					range-thickness[0]-thickness[1]
				);
				double d2d3ssd_energy = calculator.Energy(
					range-thickness[0]
				);
				double d3_energy = d3ssd_energy - ssd_energy;
				double d2_energy = d2d3ssd_energy - d3ssd_energy;
				double d1_energy = full_energy - d2d3ssd_energy;
				d1d2_pass_de_list.push_back(d1_energy);
				d1d2_pass_e_list.push_back(d2_energy);
				d2d3_pass_de_list.push_back(d2_energy);
				d2d3_pass_e_list.push_back(d3_energy);
			}
		}

		size_t points = d1d2_stop_de_list.size();
		double valbeg = d1d2_stop_de_list[1] - d1d2_stop_de_list[0];
		valbeg /= d1d2_stop_e_list[1] - d1d2_stop_e_list[0];
		double valend = d1d2_stop_de_list[points-1] - d1d2_stop_de_list[points-2];
		valend /= d1d2_stop_e_list[points-1] - d1d2_stop_e_list[points-2];

		// dE-E function, get dE from E
		TSpline3 *d1d2_stop_de_e = new TSpline3(
			"de_e_stop_0",
			d1d2_stop_e_list.data(), d1d2_stop_de_list.data(), points,
			"b1e1",
			valbeg, valend
		);

		// reverse the list, TSpline3 requests increasing x position
		std::reverse(d1d2_stop_e_list.begin(), d1d2_stop_e_list.end());
		std::reverse(d1d2_stop_de_list.begin(), d1d2_stop_de_list.end());
		// calculate valbeg and valend
		double tmp = 1.0 / valbeg;
		valbeg = 1.0 / valend;
		valend = tmp;
		// E-dE funcion, get E from dE
		TSpline3 *d1d2_stop_e_de = new TSpline3(
			"e_de_stop_0",
			d1d2_stop_de_list.data(), d1d2_stop_e_list.data(), points,
			"b1e1",
			valbeg, valend
		);

		// reverse the list, TSpline3 requests increasing x position
		std::reverse(d1d2_pass_e_list.begin(), d1d2_pass_e_list.end());
		std::reverse(d1d2_pass_de_list.begin(), d1d2_pass_de_list.end());
		points = d1d2_pass_de_list.size();
		valbeg = d1d2_pass_de_list[1] - d1d2_pass_de_list[0];
		valbeg /= d1d2_pass_e_list[1] - d1d2_pass_e_list[0];
		valend = d1d2_pass_de_list[points-1] - d1d2_pass_de_list[points-2];
		valend /= d1d2_pass_e_list[points-1] - d1d2_pass_e_list[points-2];

		// dE-E function, get dE from E
		TSpline3 *d1d2_pass_de_e = new TSpline3(
			"de_e_pass_0",
			d1d2_pass_e_list.data(), d1d2_pass_de_list.data(), points,
			"b1e1",
			valbeg, valend
		);

		// calculate valbeg and valend
		tmp = 1.0 / valbeg;
		valbeg = 1.0 / valend;
		valend = tmp;
		// E-dE funcion, get E from dE
		TSpline3 *d1d2_pass_e_de = new TSpline3(
			"e_de_pass_0",
			d1d2_pass_de_list.data(), d1d2_pass_e_list.data(), points,
			"b1e1",
			valbeg, valend
		);


		points = d2d3_stop_de_list.size();
		valbeg = d2d3_stop_de_list[1] - d2d3_stop_de_list[0];
		valbeg /= d2d3_stop_e_list[1] - d2d3_stop_e_list[0];
		valend = d2d3_stop_de_list[points-1] - d2d3_stop_de_list[points-2];
		valend /= d2d3_stop_e_list[points-1] - d2d3_stop_e_list[points-2];

		// dE-E function, get dE from E
		TSpline3 *d2d3_stop_de_e = new TSpline3(
			"de_e_stop_1",
			d2d3_stop_e_list.data(), d2d3_stop_de_list.data(), points,
			"b1e1",
			valbeg, valend
		);

		// reverse the list, TSpline3 requests increasing x position
		std::reverse(d2d3_stop_e_list.begin(), d2d3_stop_e_list.end());
		std::reverse(d2d3_stop_de_list.begin(), d2d3_stop_de_list.end());
		// calculate valbeg and valend
		tmp = 1.0 / valbeg;
		valbeg = 1.0 / valend;
		valend = tmp;
		// E-dE funcion, get E from dE
		TSpline3 *d2d3_stop_e_de = new TSpline3(
			"e_de_stop_1",
			d2d3_stop_de_list.data(), d2d3_stop_e_list.data(), points,
			"b1e1",
			valbeg, valend
		);


		// reverse the list, TSpline3 requests increasing x position
		std::reverse(d2d3_pass_e_list.begin(), d2d3_pass_e_list.end());
		std::reverse(d2d3_pass_de_list.begin(), d2d3_pass_de_list.end());
		points = d2d3_pass_de_list.size();
		valbeg = d2d3_pass_de_list[1] - d2d3_pass_de_list[0];
		valbeg /= d2d3_pass_e_list[1] - d2d3_pass_e_list[0];
		valend = d2d3_pass_de_list[points-1] - d2d3_pass_de_list[points-2];
		valend /= d2d3_pass_e_list[points-1] - d2d3_pass_e_list[points-2];

		// dE-E function, get dE from E
		TSpline3 *d2d3_pass_de_e = new TSpline3(
			"de_e_pass_1",
			d2d3_pass_e_list.data(), d2d3_pass_de_list.data(), points,
			"b1e1",
			valbeg, valend
		);

		// calculate valbeg and valend
		tmp = 1.0 / valbeg;
		valbeg = 1.0 / valend;
		valend = tmp;
		// E-dE funcion, get E from dE
		TSpline3 *d2d3_pass_e_de = new TSpline3(
			"e_de_pass_1",
			d2d3_pass_de_list.data(), d2d3_pass_e_list.data(), points,
			"b1e1",
			valbeg, valend
		);

		output_file.cd();
		d1d2_stop_de_e->Write("de_e_stop_0");
		d1d2_stop_e_de->Write("e_de_stop_0");
		d1d2_pass_de_e->Write("de_e_pass_0");
		d1d2_pass_e_de->Write("e_de_pass_0");
		d2d3_stop_de_e->Write("de_e_stop_1");
		d2d3_stop_e_de->Write("e_de_stop_1");
		d2d3_pass_de_e->Write("de_e_pass_1");
		d2d3_pass_e_de->Write("e_de_pass_1");
		output_file.Close();
	}
	return 0;
}


double D2EnergyCalculator::Energy(
	unsigned short layer,
	double delta_energy,
	double theta,
	bool stop
) const {
	// epsilon, stop if reach this standard
	double eps = 1e-4;
	// thickness of detector
	double thick1 = t0_thickness[layer];
	double thick2 = t0_thickness[layer+1];
	// // sensitive thickness in Si
	// double si_sensitive_thick = thick - dead_si_thick * 2;
	// cos(theta)
	double cos_theta = cos(theta);

	// lower bound of the total energy
	// double lower_bound = stop
	// 	? si_calculator_.Energy(thick1 / cos_theta)
	// 	: si_calculator_.Energy((thick1+thick2) / cos_theta);
	double lower_bound = si_calculator_.Energy(thick1 / cos_theta);
	// something went wrong if the Si energy is larger than the lower bound
	if (delta_energy > lower_bound) return -1e5;

	// check stop
	double stop_delta_energy =
		si_calculator_.Energy((thick1+thick2)/cos_theta)
		- si_calculator_.Energy(thick2 / cos_theta);
	if (delta_energy < stop_delta_energy) stop = false;

	// uppper bound of the total energy
	double upper_bound = stop
		? si_calculator_.Energy((thick1+thick2) / cos_theta)
		: 500.0;

	// max iteration to avoid infinite loop
	const int max_iteration = 100;
	int current_iteration = 0;
	// current energy loss in sensitive layer
	double current_si_energy = 1000.0;
	// binary search for the total energy between the lower and upper bound
	while (fabs(current_si_energy - delta_energy) > eps) {
		double current_total_energy = (upper_bound + lower_bound) / 2.0;
		// incident range in Si under current total energy
		double si_range = si_calculator_.Range(current_total_energy);
		// residual range after the sensitive layer
		si_range -= thick1 / cos_theta;
		// residual energy after sensitive layer
		double residual_energy = si_calculator_.Energy(si_range);
		// loss energy under current total energy
		current_si_energy = current_total_energy - residual_energy;
		// change the upper or lower bound
		if (current_si_energy > delta_energy) {
			// current Si energy loss too large, increase lower bound
			lower_bound = current_total_energy;
		} else {
			// current Si energy loss too large, decrease upper bound
			upper_bound = current_total_energy;
		}
		// avoid infinite loop
		++current_iteration;
		if (current_iteration == max_iteration) {
			std::cout << "D2-Energy layer " << layer
				<< ", dE " << delta_energy
				<< ", theta " << theta
				<< ", stop " << (stop ? "Y" : "N") << "\n"
				<< "lower bound " << lower_bound
				<< ", upper bound " << upper_bound << "\n"  
				<< "  iteration " << current_iteration
				<< ", full E " << current_total_energy
				<< ", dE " << current_si_energy << "\n";
		}
		if (current_iteration >= max_iteration) return -2e5;
	}
	
	// get the total energy
	double energy = (upper_bound + lower_bound) / 2.0;
	// incident range in Si under total energy
	double range = si_calculator_.Range(energy);
	// residual range after the first layer
	range -= thick1 / cos_theta;
	// residual energy after the first layer
	energy = si_calculator_.Energy(range);

	range -= thick2 / cos_theta;
	if (range > 0.0) {
		// // pass the second layer, lost part of energy
		// // residual range after the second layer
		// range -= thick2 / cos_theta;
		// // there is bug in code
		// if (range < 0.0) return -3e5;
		// residual energy after the second layer
		energy -= si_calculator_.Energy(range);
	}

	return energy;
}


double D2EnergyCalculator::DeltaEnergy(
	unsigned short layer,
	double energy,
	double theta,
	bool stop
) const {
	if (stop) {
		// incident range in the second layer
		double range = si_calculator_.Range(energy);
		// total range
		range += t0_thickness[layer] / cos(theta);
		// total energy
		double result = si_calculator_.Energy(range);
		result -= energy;
		return result;
	} else {
		// first layer thickness
		double thick1 = t0_thickness[layer];
		// second layer thickness
		double thick2 = t0_thickness[layer+1];
		// cos(theta)
		double cos_theta = cos(theta);
		// lower bound of total energy
		double lower_bound = si_calculator_.Energy((thick1+thick2) / cos_theta);
		// uppper bound of total energy
		double upper_bound = 400.0;
		if (energy > lower_bound) return -1e5;

		// max iteration to avoid infinite loop
		const int max_iteration = 1000;
		int current_iteration = 0;
		// current energy loss in first layer
		double current_energy = 1000.0;
		// epsilon, stop if reach this standard
		double eps = 1e-4;
		// binary search for the total energy between the lower and upper bound
		while (fabs(current_energy - energy) > eps) {
			double current_total_energy = (upper_bound + lower_bound) / 2.0;
			// incident range in Si under current total energy
			double si_range = si_calculator_.Range(current_total_energy);
			// residual range after the first layer
			si_range -= thick1 / cos_theta;
			// residual energy after the first layer
			double residual_energy = si_calculator_.Energy(si_range);
			// residual range after the second layer
			si_range -= thick2 / cos_theta; 
			// energy after the second layer
			double pass_energy = si_calculator_.Energy(si_range);
			// lost energy in the second layer
			current_energy = residual_energy - pass_energy;
			// change the upper or lower bound
			if (current_energy > energy) {
				// current Si energy loss too large, increase lower bound
				lower_bound = current_total_energy;
			} else {
				// current Si energy loss too large, decrease upper bound
				upper_bound = current_total_energy;
			}
			// avoid infinite loop
			++current_iteration;
			if (current_iteration >= max_iteration) return -2e5;
		}

		// get the total energy
		double total_energy = (upper_bound + lower_bound) / 2.0;
		// incident range in Si under total energy
		double range = si_calculator_.Range(total_energy);
		// residual range after the first layer
		range -= thick1 / cos_theta;
		// residual energy after the first layer
		double result = total_energy - si_calculator_.Energy(range);

		// if (result > 1000) {
		// 	std::cout << "D2 DeltaEnergy, layer " << layer
		// 		<< ", energy " << energy
		// 		<< ", theta " << theta
		// 		<< ", stop " << (stop ? "Y" : "N") << "\n"
		// 		<< "  lower bound " << lower_bound
		// 		<< ", upper bound " << upper_bound
		// 		<< ", iteration " << current_iteration
		// 		<< ", current energy " << current_energy
		// 		<< "\n"; 
		// }

		return result;
	}
}


}	// elc

}	// ribll