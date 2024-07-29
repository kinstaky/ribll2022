#include <cmath>
#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <Math/Vector3D.h>

#include "include/event/generate_event.h"
#include "include/event/threebody_info_event.h"

using namespace ribll;

int main() {
	// input file name
	TString input_file_name = TString::Format(
		"%s%sgenerate-0002.root", kGenerateDataPath, kSimulateDir
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
	GenerateEvent generate;
	// setup input branches
	generate.SetupInput(ipt);

	// 3body information file name
	TString info_file_name = TString::Format(
		"%s%sthreebody-sim-0002.root", kGenerateDataPath, kInformationDir
	);
	// add friend
	ipt->AddFriend("info=tree", info_file_name);
	// information event
	ThreeBodyInfoEvent info;
	// setup input branches
	info.SetupInput(ipt, "info.");

	// output file name
	TString output_file_name = TString::Format(
		"%s%scheck-momentum-0002.root", kGenerateDataPath, kSimulateDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// histogram of residual momentum's length
	TH1F hist_residual_p("hrp", "residual momentum value", 1000, -1, 1);
	// output tree
	TTree opt("tree", "check momentum");
	// output data
	double residual_pr, residual_px, residual_py, residual_pz;
	double h2_kinetic_difference;
	// setup output branches
	opt.Branch("residual_pr", &residual_pr, "rpr/D");
	opt.Branch("residual_px", &residual_px, "rpx/D");
	opt.Branch("residual_py", &residual_py, "rpy/D");
	opt.Branch("residual_pz", &residual_pz, "rpz/D");
	opt.Branch("h2_kinetic_difference", &h2_kinetic_difference, "h2dk/D");

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries 
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Start to check momentum   0%%");
	fflush(stdout);
	// loop
	for (long long entry = 0; entry < entries; ++entry) {
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get data
		ipt->GetEntry(entry);

		// 14C direction
		ROOT::Math::XYZVector c14_direction_gen(
			sin(generate.beam_theta) * cos(generate.beam_phi),
			sin(generate.beam_theta) * sin(generate.beam_phi),
			cos(generate.beam_theta)
		);
		// 14C momentum value
		double c14_momentum_gen = MomentumFromKinetic(
			mass_14c, generate.beam_kinetic_in_target
		);
		// 14C momentum vector
		ROOT::Math::XYZVector c14_p_gen =
			c14_direction_gen * c14_momentum_gen;

		// 10Be direction
		ROOT::Math::XYZVector be10_direction_gen(
			generate.fragment_x[0] - generate.target_x,
			generate.fragment_y[0] - generate.target_y,
			generate.fragment_z[0]
		);
		be10_direction_gen = be10_direction_gen.Unit();
		// 10Be momentum value
		double be10_momentum_gen = MomentumFromKinetic(
			mass_10be, generate.fragment_kinetic_in_target[0]
		);
		// 10Be momentum vector
		ROOT::Math::XYZVector be10_p_gen =
			be10_direction_gen * be10_momentum_gen;

		// 4He direction
		ROOT::Math::XYZVector he4_direction_gen(
			generate.fragment_x[1] - generate.target_x,
			generate.fragment_y[1] - generate.target_y,
			generate.fragment_z[1]
		);
		he4_direction_gen = he4_direction_gen.Unit();
		// 4He momentum value
		double he4_momentum_gen = MomentumFromKinetic(
			mass_4he, generate.fragment_kinetic_in_target[1]
		);
		// 4He momentum vector
		ROOT::Math::XYZVector he4_p_gen =
			he4_direction_gen * he4_momentum_gen;

		// 2H direction
		ROOT::Math::XYZVector h2_direction_gen(
			generate.rx - generate.target_x,
			generate.ry - generate.target_y,
			generate.rz
		);
		h2_direction_gen = h2_direction_gen.Unit();
		// 2H momentum value
		double h2_momentum_gen = MomentumFromKinetic(
			mass_2h, generate.recoil_kinetic_in_target
		);
		// 2H momentum vector
		ROOT::Math::XYZVector h2_p_gen =
			h2_direction_gen * h2_momentum_gen;

		// residual momentum
		ROOT::Math::XYZVector residual_p_gen =
			h2_p_gen + he4_p_gen + be10_p_gen - c14_p_gen;

		hist_residual_p.Fill(residual_p_gen.R());


		residual_pr = residual_p_gen.R();
		residual_px = residual_p_gen.X();
		residual_py = residual_p_gen.Y();
		residual_pz = residual_p_gen.Z();
	 
		// 10Be direction
		ROOT::Math::XYZVector be10_direction_info(
			info.be_x[0] - generate.target_x,
			info.be_y[0] - generate.target_y,
			100.0
		);
		be10_direction_info = be10_direction_info.Unit();
		// 10Be momentum value
		double be10_momentum_info = MomentumFromKinetic(
			mass_10be, info.t0_energy[0]
		);
		// 10Be momentum vector
		ROOT::Math::XYZVector be10_p_info =
			be10_direction_info * be10_momentum_info;

		// 4He direction
		ROOT::Math::XYZVector he4_direction_info(
			info.he_x[0] - generate.target_x,
			info.he_y[0] - generate.target_y,
			100.0
		);
		he4_direction_info = he4_direction_info.Unit();
		// 4He momentum value
		double he4_momentum_info = MomentumFromKinetic(
			mass_4he, info.t0_energy[1]
		);
		// 4He momentum vector
		ROOT::Math::XYZVector he4_p_info =
			he4_direction_info * he4_momentum_info;

		// calculate 2H momentum
		ROOT::Math::XYZVector h2_p_calc = c14_p_gen - be10_p_info - he4_p_info;
		// calculate 2H kinetic
		double h2_kinetic_calc =
			sqrt(pow(h2_p_calc.R(), 2.0) + pow(mass_2h, 2.0)) - mass_2h;

		h2_kinetic_difference =
			h2_kinetic_calc - generate.recoil_kinetic_after_target;

		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");


	// save histograms
	hist_residual_p.Write();
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}