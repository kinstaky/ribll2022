#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TString.h>
#include <Math/Vector3D.h>
#include <TRandom3.h>

#include "include/event/generate_event.h"
#include "include/calculator/target_energy_calculator.h"
#include "include/ppac_track.h"

using namespace ribll;
using ROOT::Math::XYZVector;

int main(int argc, char **argv) {
	// parse arguments
	int run = 2;
	if (argc > 1) {
		run = atoi(argv[1]);
	}

	// input file name
	TString input_file_name = TString::Format(
		"%s%sgenerate-%04d.root",
		kGenerateDataPath,
		kSimulateDir,
		run
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "[Error] Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// input event
	GenerateEvent event;
	event.SetupInput(ipt);


	// output file name
	TString output_file_name = TString::Format(
		"%s%stwobody-%04d.root",
		kGenerateDataPath,
		kSimulateDir,
		run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// Q value spectrum
	TH1F hist_q("hq", "two body Q", 300, -23, -8);
	// output tree
	TTree opt("tree", "two body spectrum");
	// output data
	int valid;
	double tx, ty;
	double bek, hek, dk, ck;
	double q;
	// setup output branches
	opt.Branch("valid", &valid, "v/I");
	opt.Branch("target_x", &tx, "tx/D");
	opt.Branch("target_y", &ty, "ty/D");
	opt.Branch("be_kinetic", &bek, "bek/D");
	opt.Branch("he_kinetic", &hek, "hek/D");
	opt.Branch("d_kinetic", &dk, "dk/D");
	opt.Branch("c_kinetic", &ck, "ck/D");
	opt.Branch("q", &q, "q/D");


	// elc::TargetEnergyCalculator c14_target("14C", "CD2", 9.53);
	// elc::TargetEnergyCalculator be10_target("10Be", "CD2", 9.53);
	// elc::TargetEnergyCalculator he4_target("4He", "CD2", 9.53);
	// elc::TargetEnergyCalculator h2_target("2H", "CD2", 9.53);

	TRandom3 generator(2853928);

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Rebuilding two body spectrum   0%%");
	fflush(stdout);
	// loop
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get data
		ipt->GetEntry(entry);
		// initialize
		valid = 0;

		tx = event.target_x + generator.Gaus(0.0, 0.65);
		ty = event.target_y + generator.Gaus(0.0, 0.65);

		// get Be information
		// kinetic energy
		bek = event.fragment_kinetic_in_target[0];
		// bek = be10_target.Energy(
		// 	-0.5, bek
		// );
		// direction
		XYZVector be_direction(
			event.fragment_x[0] - tx,
			event.fragment_y[0] - ty,
            event.fragment_z[0]
		);
		be_direction = be_direction.Unit();
		// momentum value
		double be_momentum = MomentumFromKinetic(mass_10be, bek);
		// momentum vector
		XYZVector be_p = be_direction * be_momentum;

        // get He information
		//kinetic energy
		hek = event.fragment_kinetic_in_target[1];
		// hek = he4_target.Energy(
        //     -0.5, hek
        // );
		// direction
		XYZVector he_direction(
            event.fragment_x[1] - tx,
            event.fragment_y[1] - ty,
            event.fragment_z[1]
        );
        he_direction = he_direction.Unit();
		// calculate momentum value
        double he_momentum = MomentumFromKinetic(mass_4he, hek);
        // momentum
        XYZVector he_p = he_direction * he_momentum;

		// get C information
		// kinetic energy
		ck = event.beam_kinetic_in_target;
		// ck = c14_target.Energy(
		// 	0.5, ck
		// );
		// momentum value
		double c_momentum = MomentumFromKinetic(mass_14c, ck);
		// direction
		XYZVector c_direction(
			sin(event.beam_theta) * cos(event.beam_phi),
			sin(event.beam_theta) * sin(event.beam_phi),
			cos(event.beam_theta)
		);
		c_direction = c_direction.Unit();
		// momentum vector
		XYZVector c_p = c_direction * c_momentum;

		// calculate recoil deutron
		// monmentum vector
		XYZVector d_p = c_p - be_p - he_p;
		// kinetic energy
		dk = KineticFromMomentum(run == 3 ? mass_1h : mass_2h, d_p.R());

		// Q value
		q = bek + hek + dk - ck;

        // fill histogram
        if (valid == 0) hist_q.Fill(q);
        // fill tree
        opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save
	opf.cd();
	hist_q.Write();
	opt.Write();
	// close file
	opf.Close();
	ipf.Close();
	return 0;
}
