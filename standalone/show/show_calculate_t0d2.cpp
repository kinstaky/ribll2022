#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "include/event/threebody_info_event.h"
#include "include/calculator/d2_energy_calculator.h"

using namespace ribll;

int main() {
	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody.root", kGenerateDataPath, kInformationDir
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// input event
	ThreeBodyInfoEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// output file name
	TString output_file_name = TString::Format(
		"%s%scalculate-t0d2.root", kGenerateDataPath, kShowDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "calculate T0D2 energy");
	// output data
	double d1_energy[2];
	double d2_energy[2];
	double d3_energy[2];
	double calc_d2_from_d1[2];
	double calc_d2_from_d3[2];
	// setup output branches
	opt.Branch("d1_energy", d1_energy, "d1e[2]/D");
	opt.Branch("d2_energy", d2_energy, "d2e[2]/D");
	opt.Branch("d3_energy", d3_energy, "d3e[2]/D");
	opt.Branch("calc_d2_from_d1", calc_d2_from_d1, "cd2ef1[2]/D");
	opt.Branch("calc_d2_from_d3", calc_d2_from_d3, "cd2ef3[2]/D");

	elc::D2EnergyCalculator be_calculator("10Be");
	elc::D2EnergyCalculator he_calculator("4He");

	// loop to calculate T0D2 energy
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);
		d1_energy[0] = t0_param[0][0] + t0_param[0][1]*event.be_channel[0];
		d1_energy[1] = t0_param[0][0] + t0_param[0][1]*event.he_channel[0];
		d2_energy[0] = t0_param[1][0] + t0_param[1][1]*event.be_channel[1];
		d2_energy[1] = t0_param[1][0] + t0_param[1][1]*event.he_channel[1];
		d3_energy[0] = t0_param[2][0] + t0_param[2][1]*event.be_channel[2];
		d3_energy[1] = t0_param[2][0] + t0_param[2][1]*event.he_channel[2];
		// if (event.layer[0] == 1) {
		// 	calc_d2_from_d1[0] = be_calculator.Energy(
		// 		0, d1_energy[0], event.be_theta, true
		// 	);
		// 	calc_d2_from_d3[0] = -1e5;
		// } else if (event.layer[0] >= 2) {
		// 	calc_d2_from_d1[0] = be_calculator.Energy(
		// 		0, d1_energy[0], event.be_theta, false
		// 	);
		// 	calc_d2_from_d3[0] = be_calculator.DeltaEnergy(
		// 		1, d3_energy[0], event.be_theta, event.layer[0]==2
		// 	);
		// }
		// if (event.layer[1] == 1) {
		// 	calc_d2_from_d1[1] = he_calculator.Energy(
		// 		0, d1_energy[1], event.he_theta, true
		// 	);
		// 	calc_d2_from_d3[1] = -1e5;
		// } else if (event.layer[1] >= 2) {
		// 	calc_d2_from_d1[1] = he_calculator.Energy(
		// 		0, d1_energy[1], event.he_theta, false
		// 	);
		// 	calc_d2_from_d3[1] = he_calculator.DeltaEnergy(
		// 		1, d3_energy[1], event.he_theta, event.layer[1]==2
		// 	);
		// }
		opt.Fill();
	}

	// save tree and close file
	opf.cd();
	opt.Write();
	opf.Close();
	// close file
	ipf.Close();
	return 0;
}

