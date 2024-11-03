#include <iostream>

#include <TH1F.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TChain.h>
#include <Math/Vector3D.h>

#include "include/event/particle_event.h"
#include "include/event/ta_event.h"

using namespace ribll;

int main() {
	// T0 chain
	TChain t0_chain("t0", "t0");
	// TAF chain
	TChain taf_chain[6] = {
		TChain{"taf0", "taf0"},
		TChain{"taf1", "taf1"},
		TChain{"taf2", "taf2"},
		TChain{"taf3", "taf3"},
		TChain{"taf4", "taf4"},
		TChain{"taf5", "taf5"}
	};
	// TAF particle chain
	TChain taf_particle_chain[6] = {
		TChain{"taf0p", "taf0p"},
		TChain{"taf1p", "taf1p"},
		TChain{"taf2p", "taf2p"},
		TChain{"taf3p", "taf3p"},
		TChain{"taf4p", "taf4p"},
		TChain{"taf5p", "taf5p"},
	};

	for (int run = 618; run <= 716; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		if (run > 716 && run < 739) continue;
		// T0 chain
		t0_chain.AddFile(TString::Format(
			"%s%st0-particle-ta-%04d.root/tree",
			kGenerateDataPath,
			kParticleDir,
			run
		));
		// TAF chain
		for (int i = 0; i < 6; ++i) {
			taf_chain[i].AddFile(TString::Format(
				"%s%staf%d-telescope-ta-%04d.root/tree",
				kGenerateDataPath,
				kTelescopeDir,
				i,
				run
			));
		}
		// TAF particle chain
		for (int i = 0; i < 6; ++i) {
			taf_particle_chain[i].AddFile(TString::Format(
				"%s%staf%d-particle-ta-v2-%04d.root/tree",
				kGenerateDataPath,
				kParticleDir,
				i,
				run
			));
		}
	}
	// add friends
	for (int i = 0; i < 6; ++i) {
		t0_chain.AddFriend(taf_chain+i, TString::Format("taf%d", i));
		t0_chain.AddFriend(taf_particle_chain+i, TString::Format("taf%dp", i));
	}
	// input T0 events
	ParticleEvent t0_event;
	// input TAF events
	TaEvent taf_event[6];
	// input TAF particle events
	ParticleEvent taf_particle_event[6];
	// setup input branches
	t0_event.SetupInput(&t0_chain);
	for (int i = 0; i < 6; ++i) {
		taf_event[i].SetupInput(
			&t0_chain, TString::Format("taf%d.", i).Data()
		);
		taf_particle_event[i].SetupInput(
			&t0_chain, TString::Format("taf%dp.", i).Data()
		);
	}

	// output file name
	TString output_file_name = TString::Format(
		"%s%staf-coincident-number.root",
		kGenerateDataPath,
		kShowDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "TAF coincident number");
	// output data
	int front_strip, mass, taf_index;
	opt.Branch("taf_index", &taf_index, "tafi/I");
	opt.Branch("front_strip", &front_strip, "fs/I");
	opt.Branch("mass", &mass, "A/I");

	int conflict_number = 0;

	// total number of entries
	long long entries = t0_chain.GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling TAF pid   0%%");
	fflush(stdout);
	// loop to fill pid histogram
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		t0_chain.GetEntry(entry);

		// check Be and He
		int he_index = -1;
		int be_index = -1;
		for (int i = 0; i < t0_event.num; ++i) {
			if (t0_event.charge[i] == 2 && t0_event.mass[i] == 4) he_index = i;
			if (t0_event.charge[i] == 4 && t0_event.mass[i] == 10) be_index = i;
		}
		if (he_index < 0 || be_index < 0) continue;
		// check possible deutron
		int taf_number = 0;
		taf_index = -1;
		for (int i = 0; i < 6; ++i) {
			if (taf_particle_event[i].num == 0) continue;
			if (taf_event[i].num == 0) continue;
			taf_index = i;
			++taf_number;
		}
		if (taf_index < 0) continue;
		if (taf_number > 1) {
			++conflict_number;
			continue;
		}
		front_strip = taf_event[taf_index].front_strip[0];
		mass = taf_particle_event[taf_index].mass[0];

		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	std::cout << "Multiple deutron " << conflict_number << " / "
		<< opt.GetEntries() << " "
		<< double(conflict_number) / opt.GetEntries() << "\n";


	// save
	opf.cd();
	opt.Write();
	opf.Close();

	return 0;
}