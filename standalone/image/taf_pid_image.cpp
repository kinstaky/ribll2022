#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>

#include "include/event/dssd_event.h"
#include "include/event/ta_event.h"

using namespace ribll;

int main() {
	// T0D1 chain
	TChain d1_chain("t0d1", "t0d1");
	// TAF chain
	TChain taf_chain[6] = {
		TChain{"taf0", "taf0"},
		TChain{"taf1", "taf1"},
		TChain{"taf2", "taf2"},
		TChain{"taf3", "taf3"},
		TChain{"taf4", "taf4"},
		TChain{"taf5", "taf5"}
	};
	for (unsigned int run = 628; run <= 746; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		if (run > 716 && run < 739) continue;
		d1_chain.AddFile(TString::Format(
			"%s%st0d1-merge-ta-%04d.root/tree",
			kGenerateDataPath,
			kMergeDir,
			run
		));
		for (int i = 0; i < 6; ++i) {
			taf_chain[i].AddFile(TString::Format(
				"%s%staf%d-telescope-ta-%04d.root/tree",
				kGenerateDataPath,
				kTelescopeDir,
				i,
				run
			));
		}
	}
	// add friends
	for (int i = 0; i < 6; ++i) {
		d1_chain.AddFriend(taf_chain+i, TString::Format("taf%d", i));
	}
	// input T0D1 event
	DssdMergeEvent d1_event;
	// input TAF event
	TaEvent taf_event[6];
	// setup input branches
	d1_event.SetupInput(&d1_chain);
	for (int i = 0; i < 6; ++i) {
		taf_event[i].SetupInput(
			&d1_chain, TString::Format("taf%d.", i).Data()
		);
	}

	// output histograms
	TH2F pid_all_strips[12];
	TH2F pid_single_strip[12][16];
	for (int i = 0; i < 12; ++i) {
		pid_all_strips[i] = TH2F(
			TString::Format("hc%d", i),
			TString::Format("TAF%d-CsI%d all strips PID", i/2, i),
			1000, 0, 100, 1000, 0, 20
		);
		for (int j = 0; j < 16; ++j) {
			pid_single_strip[i][j] = TH2F(
				TString::Format("hc%ds%d", i, j),
				TString::Format("TAF%d-CsI%d strip %d PID", i/2, i, j),
				1000, 0, 100, 1000, 0, 20
			);
		}
	}

	// total number of entries
	long long entries = d1_chain.GetEntries();
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
		d1_chain.GetEntry(entry);
		if (d1_event.hit == 0) continue;

		for (int i = 0; i < 6; ++i) {
			if (taf_event[i].num == 0) continue;
			if (
				taf_event[i].flag[0] != 0x3
				&& taf_event[i].flag[0] != 0x5
			) continue;
			int csi_index = taf_event[i].flag[0] == 0x3 ? i*2 : i*2+1;
			// TAFD energy
			double de = taf_event[i].energy[0][0];
			// TAF CsI energy
			double e = taf_event[i].energy[0][1];
			// calibrate CsI energy
			double a0 = power_csi_param[csi_index][0];
			double a1 = power_csi_param[csi_index][1];
			double a2 = power_csi_param[csi_index][2];
			e = pow((e-a2)/a0, 1.0/a1);
			// fill to total histogram
			pid_all_strips[csi_index].Fill(e, de);
			// fill to single strip histogram
			unsigned short &fs = taf_event[i].front_strip[0];
			pid_single_strip[csi_index][fs].Fill(e, de);
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// canvas
	TCanvas *c1 = new TCanvas("c1", "c1", 1920, 1080);
	c1->cd();

	// output pdf file name
	TString pdf_file_name = TString::Format(
		"%s%staf-pid.pdf",
		kGenerateDataPath, kImageDir
	);
	// print to pdf
	c1->Print(pdf_file_name+"[");
	for (int i = 0; i < 12; ++i) {
		pid_all_strips[i].Draw("colz");
		c1->Print(pdf_file_name);
	}
	for (int i = 0; i < 12; ++i) {
		for (int j = 0; j < 16; ++j) {
			pid_single_strip[i][j].Draw("colz");
			c1->Print(pdf_file_name);
		}
	}
	c1->Print(pdf_file_name+"]");

	return 0;
}