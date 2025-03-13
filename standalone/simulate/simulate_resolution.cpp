#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <Math/Vector3D.h>

#include "include/event/generate_event.h"
#include "include/event/channel_v2_event.h"

using namespace ribll;

int main() {
	// run number
	int run = 2;

	// input file name
	TString channel_file_name = TString::Format(
		"%s%sC14-10Be-4He-2H-v2-sim.root",
		kGenerateDataPath,
		kChannelDir
	);
	// channel file
	TFile channel_file(channel_file_name, "read");
	// input tree
	TTree *channel_tree = (TTree*)channel_file.Get("tree");
	if (!channel_tree) {
		std::cerr << "Error: Get tree from "
			<< channel_file_name << " failed.\n";
        return -1;
	}
	// setup input branches
	ChannelV2Event channel;
    channel.SetupInput(channel_tree);

	// generate file name
	TString generate_file_name = TString::Format(
        "%s%sgenerate-%04d.root",
        kGenerateDataPath,
        kSimulateDir,
		run
    );
    // generate file
    TFile generate_file(generate_file_name, "read");
    // input tree
    TTree *generate_tree = (TTree*)generate_file.Get("tree");
	if (!generate_tree) {
		std::cerr << "Error: Get tree from "
			<< generate_file_name << " failed.\n";
		return -1;
	}
	GenerateEvent generate;
	// setup input branches
	generate.SetupInput(generate_tree);

	// output file name
	TString output_file_name = TString::Format(
		"%s%sresolution-%04d.root", kGenerateDataPath, kSimulateDir, run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// histogram of excited energy difference
	TH1F hist_resolution[3][100];
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 100; ++j) {
			hist_resolution[i][j] = TH1F(
				TString::Format("hs%di%d", i, j), "energy difference", 100, -2, 2
			);
		}
	}
	// resolution graph
	TGraph graph_resolution[3];


	// if (ipt->GetEntries() != 3'000'000) {
	// 	std::cerr << "Error: entries != 3,000,000\n";
	// 	return -1;
	// }
	// 1/100 of total entries, for showing process
	long long entry100 = 30'000;
	// show start
	printf("Getting resolution   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < channel_tree->GetEntries(); ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get data
		channel_tree->GetEntry(entry);
		if (channel.valid != 0) continue;
		generate_tree->GetEntry(channel.entry);

		// 10Be state
		int state = channel.entry / 1'000'000;
		// excited energy part
		int i = (channel.entry % 1'000'000) / 10'000;


		ROOT::Math::XYZVector be_direction = ROOT::Math::XYZVector(
			channel.fragment_x[0] - channel.tx,
			channel.fragment_y[0] - channel.ty,
			100.0
		).Unit();
		double mass_excited_10be = mass_10be + generate.fragment_excited_energy;
		ROOT::Math::XYZVector bep = be_direction * MomentumFromKinetic(
			mass_excited_10be, channel.fragment_kinetic[0]
		);

		ROOT::Math::XYZVector he_direction = ROOT::Math::XYZVector(
            channel.fragment_x[1] - channel.tx,
			channel.fragment_y[1] - channel.ty,
			100.0
        ).Unit();
		ROOT::Math::XYZVector hep = he_direction * MomentumFromKinetic(
			mass_4he, channel.fragment_kinetic[1]
		);

		ROOT::Math::XYZVector xcp = bep + hep;
		double xc_energy = channel.fragment_kinetic[0] + mass_excited_10be
			+ channel.fragment_kinetic[1] + mass_4he;
		double xc_mass = sqrt(pow(xc_energy, 2.0) - xcp.Mag2());
		double xce = xc_mass - mass_14c;

		// resolution
		hist_resolution[state][i].Fill(
			xce - generate.beam_excited_energy
		);
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fit and get resolution
	for (int i = 0; i < 3; ++i) {
		double x_base = 12.02;
		if (i == 1) x_base = 15.39;
		if (i == 2) x_base = 18.20;
		for (int j = 1; j < 100; ++j) {
			TF1 *f1 = new TF1(TString::Format("fs%di%d", i, j), "gaus", -2, 2);
			hist_resolution[i][j].Fit(f1, "RQ+");
			graph_resolution[i].AddPoint(
				x_base+0.2*j, f1->GetParameter(2)
			);
		}
	}

	// write histogram
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 100; ++j) {
			hist_resolution[i][j].Write();
		}
	}
	for (int i = 0; i < 3; ++i) {
		graph_resolution[i].Write(TString::Format("g%d", i));
	}
	// close files
	opf.Close();
	generate_file.Close();
	return 0;
}