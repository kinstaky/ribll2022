#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TChain.h>
#include <TGraph.h>
#include <TH2F.h>

#include "include/event/t0_event.h"

using namespace ribll;

int main() {
	TChain chain("tree", "chain of T0 telescope events");
	for (int run = 618; run <= 627; ++run) {
		chain.AddFile(TString::Format(
			"%s%st0-telescope-%04d.root/tree",
			kGenerateDataPath, kTelescopeDir, run
		));
	}
	// input event
	T0Event event;
	// setup input branches
	event.SetupInput(&chain);

	// output file name
	TString output_file_name = TString::Format(
		"%s%spick.root",
		kGenerateDataPath, kHoleDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "hole events");
	// output data
	int front_strip, back_strip;
	double d1_energy, d2_energy;
	int type;
	// setup output branches
	opt.Branch("front_strip", &front_strip, "fs/I");
	opt.Branch("back_strip", &back_strip, "bs/I");
	opt.Branch("d1_energy", &d1_energy, "d1e/D");
	opt.Branch("d2_energy", &d2_energy, "d2e/D");
	opt.Branch("type", &type, "type/I");

	// total number of entries
	long long entries = chain.GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100;
	//show start
	printf("Processing   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get data
		chain.GetEntry(entry);

		for (int i = 0; i < event.num; ++i) {
			if (
				(event.x[i][1]+32.03)/2 == 14
				&& (event.y[i][1]+31.86)/2 == 15
				// && event.hole[i]
				&& event.layer[i] == 1
			) {
				if (event.charge[i] == 2 && event.mass[i] == 4) {
					type = 42;
				} else if (event.charge[i] == 4 && event.mass[i] == 10) {
					type = 104;
				} else if (event.charge[i] == 4 && event.mass[i] == 9) {
					type = 94;
				} else {
					type = -1;
					continue;
				}

				front_strip = int((event.x[i][1]+32.03)/2);
				back_strip = int((event.y[i][1]+31.86)/2);
				d1_energy =
					t0_param[0][0] + t0_param[0][1] * event.energy[i][0];
				d2_energy =
					t0_param[1][0] + t0_param[1][1] * event.energy[i][1];
				opt.Fill();
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save tree
	opf.cd();
	opt.Write();
	// close files
	opf.Close();
	return 0;
}