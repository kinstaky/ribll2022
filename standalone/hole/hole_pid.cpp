#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TRandom3.h>

#include "include/defs.h"
#include "include/calculator/lost_energy_calculator.h"
#include "include/calculator/d2_energy_calculator.h"

using namespace ribll;


void GetModelGraph(
	double max_energy,
	double min_energy,
	elc::LostEnergyCalculator &calculator,
	double bad_depth_start,
	double bad_depth_end,
	TGraph &pid
) {
	std::vector<double> reverse_loss;
	double energy = max_energy;
	for (int range = 0; range < 3000; ++range) {
		if (energy < 0.0) break;
		double loss = calculator.EnergyLost(energy, 1.0);
		reverse_loss.push_back(loss);
		energy -= loss;
	}

	// reverse energy loss
	std::vector<double> energy_loss;
	for (int i = reverse_loss.size()-1; i >= 0; --i) {
		energy_loss.push_back(reverse_loss[i]);
	}
	// get D1D2 energy relation
	double total_energy = 0.0;
	double d1_energy = 0.0;
	double d2_energy = 0.0;
	double energy_point = min_energy;
	for (size_t depth = 0; depth < energy_loss.size(); ++depth) {
		total_energy += energy_loss[depth];
		d1_energy += energy_loss[depth];
		if (depth >= 1010) {
			d1_energy -= energy_loss[depth-1010];
		}
		if (depth > 1010+1504) break;
		if (total_energy > energy_point && depth >= 1010) {
			energy_point += 0.1;
			d2_energy = 0.0;
			for (size_t i = 0; i <= depth-1010; ++i) {
				if (
					(depth-1010-i < 1504 && depth-1010-i > bad_depth_end)
					|| (depth-1010-i < bad_depth_start)
				) {
					d2_energy += energy_loss[i];
				}
			}
			pid.AddPoint(d1_energy, d2_energy);
		}
	}

}


int main() {
	// input file name
	TString input_file_name = TString::Format(
		"%s%spick.root", kGenerateDataPath, kHoleDir
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	// input data
	int front_strip, back_strip;
	double d1_energy, d2_energy;
	int type;
	// setup input branches
	ipt->SetBranchAddress("front_strip", &front_strip);
	ipt->SetBranchAddress("back_strip", &back_strip);
	ipt->SetBranchAddress("d1_energy", &d1_energy);
	ipt->SetBranchAddress("d2_energy", &d2_energy);
	ipt->SetBranchAddress("type", &type);

	// output file name
	TString output_file_name = TString::Format(
		"%s%spid.root", kGenerateDataPath, kHoleDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// graph to show D1D2 PID
	TGraph g_4he_pid;
	TGraph g_10be_pid;
	TH1F model_4he_diff(
		"hmd1", "difference of model and measured 4He energy", 100, -20, 20
	);
	TH1F ideal_4he_diff(
		"hid1", "difference of ideal and measured 4He energy", 100, -20, 20
	);
	TH1F model_10be_diff(
		"hmd2", "difference of model and measured 10Be energy", 100, -20, 20
	);
	TH1F ideal_10be_diff(
		"hid2", "difference of ideal and measured 10Be energy", 100, -20, 20
	);
	TH2F pid(
		"pid", "PID in D1D2 (measured)", 1000, 0, 200, 1000, 0, 200
	);
	pid.SetMarkerStyle(3);
	TH2F model_pid(
		"mpid", "PID in D1D2 with model", 1000, 0, 200, 1000, 0, 200
	);
	model_pid.SetMarkerColor(kRed);
	model_pid.SetMarkerStyle(3);
	

	// 4He energy loss calculator
	elc::LostEnergyCalculator he4_lost_calculator("4He");
	// 10Be energy loss calculator
	elc::LostEnergyCalculator be10_lost_calculator("10Be");

	// 4He T0D2 energy calculator
	elc::D2EnergyCalculator he4_d2_calculator("4He");
	// 10Be T0D2 energy calculator
	elc::D2EnergyCalculator be10_d2_calculator("10Be");


	GetModelGraph(84.0, 50.0, he4_lost_calculator, 450, 670, g_4he_pid);
	GetModelGraph(280.0, 161.0, be10_lost_calculator, 450, 670, g_10be_pid);

	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);
		if (type == 42) {
			double model_d2_energy = g_4he_pid.Eval(d1_energy);
			model_4he_diff.Fill(model_d2_energy - d2_energy);
			double ideal_d2_energy = he4_d2_calculator.Energy(
				0, d1_energy, 0.0, true
			);
			ideal_4he_diff.Fill(ideal_d2_energy - d2_energy);
			pid.Fill(d2_energy, d1_energy);
			model_pid.Fill(model_d2_energy, d1_energy);
		} else if (type == 104 || type == 94) {
			double model_d2_energy = g_10be_pid.Eval(d1_energy);
			double ideal_d2_energy = be10_d2_calculator.Energy(
				0, d1_energy, 0.0, true
			);
			if (d2_energy > 60.0) {
				model_10be_diff.Fill(model_d2_energy - d2_energy);
				ideal_10be_diff.Fill(ideal_d2_energy - d2_energy);
			}
			pid.Fill(d2_energy, d1_energy);
			model_pid.Fill(model_d2_energy, d1_energy);
		}
	}

	opf.cd();
	g_4he_pid.Write("ghe");
	g_10be_pid.Write("gbe");
	model_4he_diff.Write();
	ideal_4he_diff.Write();
	model_10be_diff.Write();
	ideal_10be_diff.Write();
	pid.Write();
	model_pid.Write();
	opf.Close();

	return 0;
}