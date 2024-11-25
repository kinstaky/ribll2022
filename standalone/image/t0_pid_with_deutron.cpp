#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>

#include "include/event/t0_event.h"
#include "include/event/particle_event.h"

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
			"%s%st0-telescope-ta-%04d.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			run
		));
		for (int i = 0; i < 6; ++i) {
			taf_chain[i].AddFile(TString::Format(
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
	}
	// input T0 event
	T0Event t0_event;
	// input TAF event
	ParticleEvent taf_event[6];
	// setup input branches
	t0_event.SetupInput(&t0_chain);
	for (int i = 0; i < 6; ++i) {
		taf_event[i].SetupInput(
			&t0_chain, TString::Format("taf%d.", i).Data()
		);
	}

	// output histograms
	TH2F d1d2_pid("d1d2", "T0 D1-D2 PID", 1000, 0, 40000, 1000, 0, 60000);
	TH2F d2d3_pid("d2d3", "T0 D2-D3 PID", 1000, 0, 30000, 1000, 0, 40000);

	// statistics
	int multiple_deutron = 0;

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

		int deutron_num = 0;
		for (int i = 0; i < 6; ++i) {
			for (int j = 0; j < taf_event[i].num; ++j) {
				if (taf_event[i].charge[i] == 1 && taf_event[i].mass[i] == 2) {
					++deutron_num;
				}	
			}
		}
		// check deutron number
		if (deutron_num > 1) {
			++multiple_deutron;
			continue;
		}
		if (deutron_num == 0) continue;
		// fill T0 particles
		for (int i = 0; i < t0_event.num; ++i) {
			if ((t0_event.flag[i] & 0x3) == 0x3) {
				d1d2_pid.Fill(t0_event.energy[i][1], t0_event.energy[i][0]);
			}
			if ((t0_event.flag[i] & 0x7) == 0x7) {
				d2d3_pid.Fill(t0_event.energy[i][2], t0_event.energy[i][1]);
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	std::cout << "Multiple deutron events " << multiple_deutron << "\n";

	// canvas
	TCanvas *c1 = new TCanvas("c1", "c1", 1920, 1080);
	c1->cd();

	// d1d2 image name
	TString d1d2_image_name = TString::Format(
		"%s%st0-pid-with-deutron-d1d2-%s.png",
		kGenerateDataPath,
		kImageDir,
		target.c_str()
	);
	d1d2_pid.Draw("colz");
	// print to pdf
	c1->Print(d1d2_image_name);

	// d2d3 image name
	TString d2d3_image_name = TString::Format(
		"%s%st0-pid-with-deutron-d2d3-%s.png",
		kGenerateDataPath,
		kImageDir,
		target.c_str()
	);
	c1->cd();
	d2d3_pid.Draw("colz");
	// print to pdf
	c1->Print(d2d3_image_name);

	return 0;
}