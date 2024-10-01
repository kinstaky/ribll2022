#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TChain.h>

#include "include/event/ta_event.h"
#include "include/event/particle_event.h"

using namespace ribll;

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " taf-index\n";
		return -1;
	}

	// TAF index
	int taf_index = atoi(argv[1]);

	// TAF chain and particle chain
	TChain taf_chain("tree", "taf chain");
	TChain particle_chain("tree", "particle chain");

	// add files
	for (int run = 618; run <= 746; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		if (run > 716 && run < 739) continue;
		taf_chain.AddFile(TString::Format(
			"%s%staf%d-telescope-ta-%04d.root",
			kGenerateDataPath,
			kTelescopeDir,
			taf_index,
			run
		));
		particle_chain.AddFile(TString::Format(
			"%s%staf%d-particle-ta-v2-%04d.root",
			kGenerateDataPath,
			kParticleDir,
			taf_index,
			run
		));
	}
	taf_chain.AddFriend(&particle_chain, "particle");

	// input data
	TaEvent taf_event;
	ParticleEvent particle_event;
	// setup input branches
	taf_event.SetupInput(&taf_chain);
	particle_event.SetupInput(&taf_chain, "particle.");

	// output file name
	TString output_file_name = TString::Format(
		"%s%staf%d-threshold.root",
		kGenerateDataPath,
		kShowDir,
		taf_index
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// threshold
	TH1F hist_tafde[16];
	TH1F hist_csie[2][16];
	TH1F hist_csie_all[2];
	for (int i = 0; i < 16; ++i) {
		hist_tafde[i] = TH1F(
			TString::Format("hd%d", i),
			TString::Format("strip %d TAFD energy", i),
			100, 0, 5
		);
		hist_csie[0][i] = TH1F(
			TString::Format("hac%d", i),
			TString::Format("strip A-%d CsI energy", i),
			100, 0, 30000
		);
		hist_csie[1][i] = TH1F(
			TString::Format("hbc%d", i),
			TString::Format("strip B-%d CsI energy", i),
			100, 0, 30000
		);
	}
	hist_csie_all[0] = TH1F("hace", "CsI-A energy", 100, 0, 30000);
	hist_csie_all[1] = TH1F("hbce", "CsI-B energy", 100, 0, 30000);

	// total number of entries
	long long entries = taf_chain.GetEntries();
	// 1/100 number of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Checking threshold   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		taf_chain.GetEntry(entry);
		if (particle_event.num == 1 && particle_event.mass[0] == 2) {
			hist_tafde[taf_event.front_strip[0]].Fill(
				taf_event.energy[0][0]
			);
			int csi_index = taf_event.flag[0] == 3 ? 0 : 1;
			hist_csie[csi_index][taf_event.front_strip[0]].Fill(
				taf_event.energy[0][1]
			);
			hist_csie_all[csi_index].Fill(taf_event.energy[0][1]);
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save histograms
	for (int i = 0; i < 16; ++i) hist_tafde[i].Write();
	for (int i = 0; i < 16; ++i) {
		hist_csie[0][i].Write();
		hist_csie[1][i].Write();
	}
	hist_csie_all[0].Write();
	hist_csie_all[1].Write();
	// close file
	opf.Close();
	return 0;
}