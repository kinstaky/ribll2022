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

/// @brief fill bragg curv into graph
/// @param[in] max_energy maximum energy
/// @param[in] eps epsilon
/// @param[in] step energy step in loop
/// @param[in] calculator lost energy calculator
/// @param[out] peak graph of bragg peak
///
void FillBraggPeak(
	double max_energy,
	elc::LostEnergyCalculator &calculator,
	TGraph &peak
) {
	// eps
	const double eps = 1e-3;
	// step
	const double step = 0.1;
	// current energy in loop
	double energy = max_energy;
	// current depth
	double depth = 0.0;
	// loop
	while (energy > eps) {
		double loss = calculator.EnergyLost(energy, step);
		if (energy < loss) {
			loss = energy;
		}
		peak.AddPoint(depth, loss);
		depth += step;
		energy -= loss;
	}
}


/// @brief a tiny event structrue for convenient
struct ModelEvent {
	double ideal_t0_energy[3];
	double t0_energy[3];
	int type;
};


/// @brief generate events and fill to output tree 
/// @param[in] energy_mean mean value of energy in generated events
/// @param[in] energy_sigma sigma value of energy in generated events 
/// @param[in] calculator lost energy calculator 
/// @param[in] type particle type, mass*10+charge, 14C:146, 10Be:104...
/// @param[out] event output event 
/// @param[out] opt output tree for storing events
/// 
void FillSimulatedEvents(
	double energy_mean,
	double energy_sigma,
	elc::LostEnergyCalculator &calculator,
	int type,
	ModelEvent &event,
	TTree &opt
) {
	TRandom3 generator(energy_mean);
	// total number of events
	const int events = 10000;
	// step
	const double step = 1.0;
	// epsilon
	const double eps = 1e-3;
	// current depth
	double depth = 0.0;
	// generate events
	for (int entry = 0; entry < events; ++entry) {
		double current_energy = generator.Gaus(energy_mean, energy_sigma);
		depth = 0.0;
		for (int i = 0; i < 3; ++i) {
			event.t0_energy[i] = 0.0;
			event.ideal_t0_energy[i] = 0.0;
		}
		event.type = type;
		while (current_energy > eps && depth < 4000.0) {
			double loss = calculator.EnergyLost(current_energy, step);
			if (current_energy < loss) {
				loss = current_energy;
			}
			if (depth < 1000.0) {
				event.t0_energy[0] += loss;
				event.ideal_t0_energy[0] += loss;
			} else if (depth < 2500.0) {
				if (depth > bad_depth_start && depth < bad_depth_end) {
					// do nothing is bad depth now
				} else {
					event.t0_energy[1] += loss;
				}
				event.ideal_t0_energy[1] += loss;
			} else {
				event.t0_energy[2] += loss;
				event.ideal_t0_energy[2] += loss;
			}
			depth += step;
			current_energy -= loss;
		}
		opt.Fill();
	}
}


void FillDepthEnergyGraph(
	double full_energy,
	elc::LostEnergyCalculator &calculator,
	TGraph &depth_energy
) {
	// energy lost in T0D1
	double d1_lost_energy = 0.0;
	// get energy loss in T0D1
	for (double range = 0.0; range < t0_thickness[0]; range += 1.0) {
		double loss = calculator.EnergyLost(full_energy-d1_lost_energy, 1.0);
		d1_lost_energy += loss;
	}

	// energy lost in T0D2
	double d2_lost_energy = 0.0;
	// energy lost in T0D2 after behind bad depth
	double d2_behind_energy = 0.0;
	// get energy lost in T0D2
	for (double range = 0.0; range < t0_thickness[1]; range += 1.0) {
		if (full_energy <= d1_lost_energy + d2_lost_energy) break;
		double loss = calculator.EnergyLost(
			full_energy-d1_lost_energy-d2_lost_energy, 1.0
		);
		d2_lost_energy += loss;
		if (range + t0_thickness[0] > bad_depth_end) {
			d2_behind_energy += loss;
		}
	}

	// current eneryg in current depth in T0D2
	double current_energy = full_energy - d1_lost_energy;
	// T0D2 meansured energy
	double energy = 0.0;
	// loop depth in T0D2
	for (
		double depth = 0.0;
		depth < bad_depth_end-t0_thickness[0];
		depth += 0.1
	) {
		double loss = calculator.EnergyLost(current_energy, 0.1);
		current_energy -= loss;
		if (current_energy < 0.0) break;
		energy += loss;
		depth_energy.AddPoint(energy+d2_behind_energy, depth);
	}
}


int main() {
	// output file name
	TString output_file_name = TString::Format(
		"%s%smodel.root",
		kGenerateDataPath, kHoleDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// graph to show Bragg peak
	TGraph g_16n_peak;
	TGraph g_14c_peak;
	TGraph g_10be_peak;
	TGraph g_9be_peak;
	TGraph g_4he_peak;
	TGraph g_14c_depth_energy;
	// output tree
	TTree opt("tree", "simulated energy lost");
	// output data
	ModelEvent event;
	// setup output branches
	opt.Branch("energy", event.t0_energy, "e[3]/D");
	opt.Branch("ideal_energy", event.ideal_t0_energy, "ie[3]/D");
	opt.Branch("type", &event.type, "type/I");

	// 16N energy loss calculator
	elc::LostEnergyCalculator n16_calculator("16N");
	// 14C energy loss calculator
	elc::LostEnergyCalculator c14_calculator("14C");
	// 10Be energy loss calculator
	elc::LostEnergyCalculator be10_calculator("10Be");
	// 9Be energy loss calculator
	elc::LostEnergyCalculator be9_calculator("9Be");
	// 4He energy loss calculator
	elc::LostEnergyCalculator he4_calculator("4He");


	// fill bragg peak for different type of particle
	FillBraggPeak(480.0, n16_calculator, g_16n_peak);
	FillBraggPeak(390.0, c14_calculator, g_14c_peak);
	FillBraggPeak(450.0, be10_calculator, g_10be_peak);
	FillBraggPeak(450.0, be9_calculator, g_9be_peak);
	FillBraggPeak(450.0, he4_calculator, g_4he_peak);

	// fill simulated events into output tree
	FillSimulatedEvents(390.0, 3.5, c14_calculator, 146, event, opt);
	FillSimulatedEvents(229.0, 31.0, be10_calculator, 104, event, opt);
	FillSimulatedEvents(210.0, 26.0, be9_calculator, 94, event, opt);
	FillSimulatedEvents(100.0, 20.0, he4_calculator, 42, event, opt);

	// fill depth energy relation graph
	FillDepthEnergyGraph(390.0, c14_calculator, g_14c_depth_energy);

	opf.cd();
	g_16n_peak.Write("gn16p");
	g_14c_peak.Write("gc14p");
	g_10be_peak.Write("gbe9p");
	g_10be_peak.Write("gbe10p");
	g_4he_peak.Write("ghe4p");
	g_14c_depth_energy.Write("gc14d");
	opt.Write();
	opf.Close();

	return 0;
}