#include <iostream>

#include <TChain.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TSpline.h>

#include "include/event/t0_event.h"
#include "include/telescope/t0.h"

using namespace ribll;

int main(int argc, char **argv) {
	if (argc < 3) {
		std::cout << "Usage: " << argv[0] << " run end_run\n"
			" run           start run number\n"
			" end_run       end run number\n";
		return -1;
	}
	int run = atoi(argv[1]);
	int end_run = atoi(argv[2]);

	// input T0 telescope chain
	TChain t0_chain("tree", "t0 tree");
	// add files
	for (int i = run; i <= end_run; ++i) {
		if (i == 628) continue;
		t0_chain.AddFile(TString::Format(
			"%s%st0-telescope-ta-%04d.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			i
		));
	}
	// input event
	T0Event t0_event;
	// setup input branches
	t0_event.SetupInput(&t0_chain);

	// output file name
	TString output_file_name = TString::Format(
		"%s%scheck-t0-pid-%04d-%04d.root",
		kGenerateDataPath,
		kCheckDir,
		run,
		end_run
	);
	// output file
	TFile output_file(output_file_name, "recreate");
	// 2D histogram of d1d2 pid
	TH2F d1d2_pid(
		"d1d2", "d1d2 #DeltaE-E pid",
		3000, 0, 300, 3000, 0, 300
	);
	// 2D histogram of d2d3 pid
	TH2F d2d3_pid(
		"d2d3", "d2d3 #DeltaE-E pid",
		2000, 0, 200, 2500, 0, 250
	);
	// 2D histogram of d3s1
	TH2F d3s1_pid(
		"d3s1", "d3s1 #DeltaE-E pid",
		1500, 0, 150, 1500, 0, 150
	);
	// 2D histogram of s1s2
	TH2F s1s2_pid(
		"s1s2", "s1s2 #DeltaE-E pid",
		2000, 0, 100, 2000, 0, 100
	);
	// 2D histogram of s2s3
	TH2F s2s3_pid(
		"s2s3", "s2s3 #DeltaE-E pid",
		2000, 0, 100, 2000, 0, 100
	);
	T0 t0(run, "ta");
	t0.ReadCalibrateParameters();

	// total number of entries
	long long entries = t0_chain.GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling T0 pid   0%%");
	fflush(stdout);
	// loop to fill pid histogram
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get entry
		t0_chain.GetEntry(entry);
		t0.CalibrateResult(t0_event);
		for (int i = 0; i < t0_event.num; ++i) {
			if ((t0_event.flag[i] & 0x3) == 0x3) {
				d1d2_pid.Fill(
					t0_event.energy[i][1], t0_event.energy[i][0]
				);
			}
			if (t0_event.flag[i] == 0x7) {
				d2d3_pid.Fill(
					t0_event.energy[i][2], t0_event.energy[i][1]
				);
				if ((t0_event.ssd_flag & 0x1) == 0x1) {
					d3s1_pid.Fill(
						t0_event.ssd_energy[0], t0_event.energy[i][2]
					);
				}
			}
			if ((t0_event.ssd_flag & 0x3) == 0x3) {
				s1s2_pid.Fill(
					t0_event.ssd_energy[1], t0_event.ssd_energy[0]
				);
			}
			if (t0_event.ssd_flag == 0x7) {
				s2s3_pid.Fill(
					t0_event.ssd_energy[2], t0_event.ssd_energy[1]
				);
			}
		}
	}
	printf("\b\b\b\b100%%\n");

	// save pid histograms
	output_file.cd();
	d1d2_pid.Write();
	d2d3_pid.Write();
	d3s1_pid.Write();
	s1s2_pid.Write();
	s2s3_pid.Write();

	std::string particles[] = {
		"3H", "3He", "4He", "6He", "6Li", "7Li", "9Be", "10Be"
	};
	for (const std::string &particle : particles) {
		TFile pid_file(
			TString::Format(
				"%s%st0-delta-%s.root",
				kGenerateDataPath,
				kEnergyCalculateDir,
				particle.c_str()
			),
			"read"
		);
		TSpline3 *curv = (TSpline3*)pid_file.Get("de_e_1");
		output_file.cd();
		curv->Write(TString::Format("pid_d2d3_%s", particle.c_str()));
		pid_file.Close();
	}
	// close files
	output_file.Close();
	return 0;
}