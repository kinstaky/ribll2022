#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>

#include "include/event/generate_event.h"

using namespace ribll;

int main() {
	// run number
	int run = 2;

	// input file name
	TString generate_file_name = TString::Format(
		"%s%sgenerate-%04d.root", kGenerateDataPath, kSimulateDir, run
	);
	// input file
	TFile generate_file(generate_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)generate_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< generate_file_name << " failed.\n";
		return -1;
	}
	// generate event
	GenerateEvent generate;
	// setup input branches
	generate.SetupInput(ipt);

	// spectrum-V2 file
	TString spectrum_file_name = TString::Format(
		"%s%sthreebody-sim-%04d-2.root", kGenerateDataPath, kSpectrumDir, run 
	);
	// add friend
	ipt->AddFriend("s=tree", spectrum_file_name);
	// PPAC flag and TAF flag and bind flag
	int ppac_flag, taf_flag, bind_flag;
	// rebuilt excited energy
	double excited_energy[4];
	// setup input branches
	ipt->SetBranchAddress("s.ppac_flag", &ppac_flag);
	ipt->SetBranchAddress("s.taf_flag", &taf_flag);
	ipt->SetBranchAddress("s.bind", &bind_flag);
	ipt->SetBranchAddress("s.excited_energy", excited_energy);

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
	
	
	if (ipt->GetEntries() != 3'000'000) {
		std::cerr << "Error: entries != 3,000,000\n";
		return -1;
	}
	// 1/100 of total entries, for showing process
	long long entry100 = 30'000;
	// show start
	printf("Getting resolution   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < 3'000'000; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get data
		ipt->GetEntry(entry);
		if ((ppac_flag & 1) == 0) continue;
		if (taf_flag) continue;
		if (bind_flag) continue;

		// 10Be state
		int state = entry / 1'000'000;
		// excited energy part
		int i = (entry % 1'000'000) / 10'000;

		// resolution
		hist_resolution[state][i].Fill(
			excited_energy[3] - generate.beam_excited_energy
		);
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fit and get resolution
	for (int i = 0; i < 3; ++i) {
		double x_base = 12.02;
		if (i == 1) x_base = 15.39;
		if (i == 2) x_base = 18.20;
		for (int j = 0; j < 100; ++j) {
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