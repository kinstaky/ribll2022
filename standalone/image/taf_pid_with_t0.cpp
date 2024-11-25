#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>

#include "include/event/particle_event.h"
#include "include/event/ta_event.h"

using namespace ribll;

bool CheckTarget(int run, const std::string target) {
	if (run == 628 || run == 723) {
		return false;
	}
	if (target == "cd2") {
		if (run >= 618 && run <= 652) return true;
		else if (run >= 675 && run <= 716) return true;
		else if (run >= 739 && run <= 746) return true;
	} else if (target == "c") {
		if (run == 656 || run == 657) return true;
		else if (run == 668 || run == 669) return true;
		else if (run >= 717 && run <= 732) return true;
	}
	return false;
}


int main(int argc, char **argv) {
	if (argc > 2 || (argc == 2 && argv[1][0] == '-' && argv[1][1] == 'h')) {
		std::cout << "Usage: " << argv[0] << " [target].\n"
			<< "  target          cd2 or c, default is cd2.\n";
		return 0;
	}

	std::string target = "cd2";
	if (argc == 2) target = std::string(argv[1]);
	if (target != "cd2" && target != "c" ) {
		std::cerr << "Error: Invalid target argument, "
			<< target << "\n";
		return -1;
	}
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
	int start_run = target == "cd2" ? 618 : 656;
	int end_run = target == "cd2" ? 746 : 732;
	for (int run = start_run; run <= end_run; ++run) {
		if (!CheckTarget(run, target)) continue;
		t0_chain.AddFile(TString::Format(
			"%s%st0-particle-ta-v2-%04d.root/tree",
			kGenerateDataPath,
			kParticleDir,
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
		t0_chain.AddFriend(taf_chain+i, TString::Format("taf%d", i));
	}
	// input T0 event
	ParticleEvent t0_event;
	// input TAF event
	TaEvent taf_event[6];
	// setup input branches
	t0_event.SetupInput(&t0_chain);
	for (int i = 0; i < 6; ++i) {
		taf_event[i].SetupInput(
			&t0_chain, TString::Format("taf%d.", i).Data()
		);
	}

	// output histograms
	TH2F pid_all_strips[12];
	// TH2F pid_single_strip[12][16];
	for (int i = 0; i < 12; ++i) {
		pid_all_strips[i] = TH2F(
			TString::Format("hc%d", i),
			TString::Format(
				"TAF%d-CsI%d all strips PID, target %s",
				i/2, i, target.c_str()
			),
			1000, 0, 100, 1000, 0, 20
		);
		// for (int j = 0; j < 16; ++j) {
		// 	pid_single_strip[i][j] = TH2F(
		// 		TString::Format("hc%ds%d", i, j),
		// 		TString::Format(
		// 			"TAF%d-CsI%d strip %d PID target %s",
		// 			i/2, i, j, target.c_str()
		// 		),
		// 		1000, 0, 100, 1000, 0, 20
		// 	);
		// }
	}

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

		bool found_be10 = false;
		bool found_he4 = false;
		for (int i = 0; i < t0_event.num; ++i) {
			if (t0_event.charge[i] == 4 && t0_event.mass[i] == 10) {
				found_be10 = true;
			} else if (t0_event.charge[i] == 2 && t0_event.mass[i] == 4) {
				found_he4 = true;
			}
		}

		if (!found_be10 || !found_he4) continue;


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
			// unsigned short &fs = taf_event[i].front_strip[0];
			// pid_single_strip[csi_index][fs].Fill(e, de);
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// canvas
	TCanvas *c1 = new TCanvas("c1", "c1", 1920, 1080);
	c1->cd();

	// output pdf file name
	TString pdf_file_name = TString::Format(
		"%s%staf-pid-with-t0-%s.pdf",
		kGenerateDataPath,
		kImageDir,
		target.c_str()
	);
	// print to pdf
	c1->Print(pdf_file_name+"[");
	for (int i = 0; i < 12; ++i) {
		pid_all_strips[i].Draw(
			target == "cd2" ? "colz" : "*"
		);
		c1->Print(pdf_file_name);
	}
	// for (int i = 0; i < 12; ++i) {
	// 	for (int j = 0; j < 16; ++j) {
	// 		pid_single_strip[i][j].Draw("colz");
	// 		c1->Print(pdf_file_name);
	// 	}
	// }
	c1->Print(pdf_file_name+"]");

	return 0;
}