#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TRandom3.h>

#include "include/defs.h"
#include "include/calculator/lost_energy_calculator.h"

using namespace ribll;

constexpr double bad_depth_start = 1400.0;
constexpr double bad_depth_end = 1680.0;


void BadDepthProcess(double &, double) {
	return;
}

int main() {

	// output file name
	TString output_file_name = TString::Format(
		"%s%shole-model.root",
		kGenerateDataPath, kShowDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// graph to show Bragg peak
	TGraph g_14c_peak;
	TGraph g_10be_peak;
	TGraph g_4he_peak;
	// output tree
	TTree opt("tree", "simulated energy lost");
	// output data
	double ideal_t0_energy[3];
	double t0_energy[3];
	int type;
	// setup output branches
	opt.Branch("energy", t0_energy, "e[3]/D");
	opt.Branch("ideal_energy", ideal_t0_energy, "ie[3]/D");
	opt.Branch("type", &type, "type/I");

	// 14C energy loss calculator
	elc::LostEnergyCalculator c14_calculator("14C");
	// 10Be energy loss calculator
	elc::LostEnergyCalculator be10_calculator("10Be");
	// 4He energy loss calculator
	elc::LostEnergyCalculator he4_calculator("4He");

	// 14C beam full energy
	const double full_energy = 390.0;

	// epsilon
	double eps = 1e-3;
	// current energy
	double energy = full_energy;
	// current depth
	double depth = 0.0;
	// step
	double step = 0.1;
	// calculate 14C 
	while (energy > eps) {
		double loss = c14_calculator.EnergyLost(energy, step);
		if (energy < loss) {
			loss = energy;
		}
		g_14c_peak.AddPoint(depth, loss);
		depth += step;
		energy -= loss;
	}

	// calculate 10Be
	energy = full_energy;
	depth = 0.0;
	while (energy > eps) {
		double loss = be10_calculator.EnergyLost(energy, step);
		if (energy < loss) {
			loss = energy;
		}
		g_10be_peak.AddPoint(depth, loss);
		depth += step;
		energy -= loss;
	}

	// calculate 4He
	energy = full_energy;
	depth = 0.0;
	while (energy > eps) {
		double loss = he4_calculator.EnergyLost(energy, step);
		if (energy < loss) {
			loss = energy;
		}
		g_4he_peak.AddPoint(depth, loss);
		depth += step;
		energy -= loss;
	}

	TRandom3 generator(2874);

	const int events = 10000;
	step = 1.0;
	// generate 14C events
	for (int i = 0; i < events; ++i) {
		double current_energy = full_energy + generator.Gaus(0.0, 3.5);
		depth = 0.0;
		for (int i = 0; i < 3; ++i) {
			t0_energy[i] = 0.0;
			ideal_t0_energy[i] = 0.0;
		}
		type = 6;
		while (current_energy > eps && depth < 4000.0) {
			double loss = c14_calculator.EnergyLost(current_energy, step);
			if (current_energy < loss) {
				loss = energy;
			}
			if (depth < 1000.0) {
				t0_energy[0] += loss;
				ideal_t0_energy[0] += loss;
			} else if (depth < 2500.0) {
				if (depth > bad_depth_start && depth < bad_depth_end) {
					BadDepthProcess(t0_energy[1], loss);
				} else {
					t0_energy[1] += loss;
				}
				ideal_t0_energy[1] += loss;
			} else {
				t0_energy[2] += loss;
				ideal_t0_energy[2] += loss;
			}
			depth += step;
			current_energy -= loss;
		}
		opt.Fill();
	}

	// generate 10Be events
	for (int i = 0; i < events; ++i) {
		double current_energy = generator.Gaus(229.0, 31.0);
		depth = 0.0;
		for (int i = 0; i < 3; ++i) {
			t0_energy[i] = 0.0;
			ideal_t0_energy[i] = 0.0;
		}
		type = 4;
		while (current_energy > eps && depth < 4000.0) {
			double loss = be10_calculator.EnergyLost(current_energy, step);
			if (current_energy < loss) {
				loss = energy;
			}
			if (depth < 1000.0) {
				t0_energy[0] += loss;
				ideal_t0_energy[0] += loss;
			} else if (depth < 2500.0) {
				if (depth > bad_depth_start && depth < bad_depth_end) {
					BadDepthProcess(t0_energy[1], loss);
				} else {
					t0_energy[1] += loss;
				}
				ideal_t0_energy[1] += loss;
			} else {
				t0_energy[2] += loss;				
				ideal_t0_energy[2] += loss;
			}
			depth += step;
			current_energy -= loss;

		}
		opt.Fill();
	}

	// generate 4He events
	for (int i = 0; i < events; ++i) {
		double current_energy = generator.Gaus(100.0, 20.0);
		depth = 0.0;
		for (int i = 0; i < 3; ++i) {
			t0_energy[i] = 0.0;
			ideal_t0_energy[i] = 0.0;
		}
		type = 2;
		while (current_energy > eps && depth < 4000.0) {
			double loss = he4_calculator.EnergyLost(current_energy, step);
			if (current_energy < loss) {
				loss = energy;
			}
			if (depth < 1000.0) {
				t0_energy[0] += loss;
				ideal_t0_energy[0] += loss;
			} else if (depth < 2500.0) {
				if (depth > bad_depth_start && depth < bad_depth_end) {
					BadDepthProcess(t0_energy[1], loss);
				} else {
					t0_energy[1] += loss;
				}
				ideal_t0_energy[1] += loss;
			} else {
				t0_energy[2] += loss;
				ideal_t0_energy[2] += loss;
			}
			depth += step;
			current_energy -= loss;
		}
		opt.Fill();
	}


	opf.cd();
	g_14c_peak.Write("gc");
	g_10be_peak.Write("gbe");
	g_4he_peak.Write("ghe");
	opt.Write();
	opf.Close();

	return 0;
}