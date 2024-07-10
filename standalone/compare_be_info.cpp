#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/threebody_info_event.h"
#include "include/calculator/d2_energy_calculator.h"

using namespace ribll;

// striaght parameters
constexpr double straight_parameters[] = {
	0.47, -0.042,	// D1D2: A, B
	0.59, -0.042	// D2D3: A, B
};

// Be mean and sigma
constexpr double be_ef_parameters[] = {
	19796.1, 250.0,		// D1D2-9Be: mean, sigma
	20757.5, 214.2,		// D1D2-10Be: mean, sigma
	23569.6, 452.53,	// D2D3-9Be: mean, sigma
	24557.7, 293.941	// D2D3-10Be: mean, sigma
};


int main() {
	// energy calculator
	elc::D2EnergyCalculator be10_calculator("10Be");

	// output file name
	TString output_file_name = TString::Format(
		"%s%scompare-be.root",
		kGenerateDataPath, kInformationDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "compare 9Be and 10Be information");
	// output data
	int mass;
	// bool valid;
	bool in_target, in_center;
	double c14_kinetic, c14_momentum;
	double c14_kinetic_sigma, c14_momentum_sigma;
	double be_ef, be_calc_ef;
	double be9_sigma, be9_calc_sigma;
	double be10_sigma, be10_calc_sigma;
	// setup output branches
	opt.Branch("mass", &mass, "A/I");
	// opt.Branch("valid", &valid, "valid/O");
	opt.Branch("in_target", &in_target, "it/O");
	opt.Branch("in_center", &in_center, "ic/O");
	opt.Branch("c14_kinetic", &c14_kinetic, "c14k/D");
	opt.Branch("c14_kinetic_sigma", &c14_kinetic_sigma, "c14ks/D");
	opt.Branch("c14_momentum", &c14_momentum, "c14p/D");
	opt.Branch("c14_momentum_sigma", &c14_momentum_sigma, "c14ps/D");
	opt.Branch("be_ef", &be_ef, "beef/D");
	opt.Branch("be_calc_ef", &be_calc_ef, "becef/D");
	opt.Branch("be9_sigma", &be9_sigma, "be9s/D");
	opt.Branch("be9_calc_sigma", &be9_calc_sigma, "be9cs/D");
	opt.Branch("be10_sigma", &be10_sigma, "be10s/D");
	opt.Branch("be10_calc_sigma", &be10_calc_sigma, "be10cs/D");



	for (mass = 9; mass <= 10; ++mass) {
		// input Be file name
		TString be_file_name = TString::Format(
			"%s%sthreebody-calc-%dBe.root",
			kGenerateDataPath, kInformationDir, mass
		);
		// input file
		TFile ipf(be_file_name, "read");
		// input tree
		TTree *ipt = (TTree*)ipf.Get("tree");
		if (!ipt) {
			std::cerr << "Error: Get tree from "
				<< be_file_name << " failed.\n";
			return -1;
		}
		// input event
		ThreeBodyInfoEvent event;
		// setup input branches
		event.SetupInput(ipt);
		
		for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
			// get data
			ipt->GetEntry(entry);

			// initialize
			// valid = false;
			in_target = false;
			in_center = false;
			c14_kinetic = 0.0;
			c14_momentum = 0.0;
			be9_sigma = 10.0;
			be10_sigma = 10.0;

			// ignore complex conditions
			if (event.taf_flag != 0) continue;
			if (event.xppac_track[0] == 0 || event.xppac_track[1] == 0) continue;
			if (event.hole[0] || event.hole[1]) continue;
			if (event.layer[0] != 1) continue;
			
			// set valid
			// valid = true;

			// check reaction position
			in_center =
				event.xptx > -10 && event.xptx < -5
				&& event.xpty > -3 && event.xpty < 3;

			// get kinetic energy and momentum
			c14_kinetic = event.c14_kinetic;
			c14_momentum = event.c14_momentum;
			if (in_center) {
				c14_kinetic_sigma = (c14_kinetic - 376.0) / 7.4;
				c14_momentum_sigma = (c14_momentum - 3147.0) / 37.3;
			} else {
				c14_kinetic_sigma = (c14_kinetic - 384.0) / 4.6;
				c14_momentum_sigma = (c14_momentum - 3184.0) / 23.7;
			}


			// get sigma in straight PID
			// 10Be direction
			ROOT::Math::XYZVector d_be(
				event.be_x[0] - event.xptx,
				event.be_y[0] - event.xpty,
				100.0
			);
			// T0D1 measured energy channel
			double d1e = event.be_channel[0];
			// T0D2 measured energy channel
			double d2e = event.be_channel[1];
			// calculated energy in T0D2
			double cd2e = be10_calculator.Energy(
				0, t0_param[0][0] + t0_param[0][1]*d1e, d_be.Theta(), true
			);
			// T0D2 calculated energy channel
			cd2e = (cd2e - t0_param[1][0]) / t0_param[1][1];
			// measured Be Ef
			be_ef =
				sqrt(d1e*d2e + straight_parameters[0]*d1e*d1e)
				+ straight_parameters[1]*d2e;
			// calculated Be Ef
			be_calc_ef =
				sqrt(d1e*cd2e + straight_parameters[0]*d1e*d1e)
				+ straight_parameters[1]*cd2e;

			// measured 9Be Ef sigma
			be9_sigma = (be_ef - be_ef_parameters[0]) / be_ef_parameters[1];
			// calculated 9Be Ef sigma
			be9_calc_sigma =
				(be_calc_ef - be_ef_parameters[0]) / be_ef_parameters[1];
			// measured 10Be Ef sigma
			be10_sigma = (be_ef - be_ef_parameters[2]) / be_ef_parameters[3];
			// calculated 10Be Ef sigma
			be10_calc_sigma =
				(be_calc_ef - be_ef_parameters[2]) / be_ef_parameters[3];


			// fill
			opt.Fill();
		}

		// close input file
		ipf.Close();
	}

	// save tree
	opf.cd();
	opt.Write();
	// close files
	opf.Close();
	
	return 0;
}