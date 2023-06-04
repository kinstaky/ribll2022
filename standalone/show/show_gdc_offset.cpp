#include <iostream>
#include <vector>

#include <TString.h>
#include <TFile.h>
#include <TGraph.h>
#include <TTree.h>

#include "include/defs.h"

using namespace ribll;

const int range = 1000;


int main(int argc, char **argv) {
	if (argc != 2) {
		std::cout << "Usage: " << argv[0] << " run" << std::endl;
		return -1;
	}

	int run = atoi(argv[1]);

	// input VME file name
	TString vme_file_name;
	vme_file_name.Form(
		"%s%s%04d.root",
		kCrate3Path, kCrate3FileName, run
	);
	// input VME file
	TFile vme_file(vme_file_name, "read");
	// input VME tree
	TTree *vme_tree = (TTree*)vme_file.Get("tree");
	if (!vme_tree) {
		std::cerr << "Error: Get tree from "
			<< vme_file_name << " failed.\n";
		return -1;
	}
	// input data
	int gmulti[2][128];
	int madc[2][32];
	// setup input branches
	vme_tree->SetBranchAddress("gmulti", gmulti);
	vme_tree->SetBranchAddress("madc", madc);
	std::vector<bool> time[32];
	std::vector<bool> energy[32];

	// total number of entries
	long long entries = vme_tree->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Reading values   0%%");
	fflush(stdout);
	// read events
	for (long long entry = 0; entry < vme_tree->GetEntriesFast(); ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		vme_tree->GetEntry(entry);

		// read values into arrays
		for (int i = 0; i < 32; ++i) {
			time[i].push_back(gmulti[1][i] > 0 ? true : false);
			size_t module_index = vtaf_front_module[i/16];
			size_t channel_index = vtaf_front_channel[i/16] + i%16;
			energy[i].push_back(
				madc[module_index][channel_index] > 200 ? true : false
			);
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// close input file
	vme_file.Close();

	// output file
	TString output_file_name;
	output_file_name.Form(
		"%s%sgdc-offset-%04d.root",
		kGenerateDataPath, kShowDir, run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output match histogram
	TGraph gm[2];

	int size = energy[0].size();
	for (int offset = 500; offset < 1000; ++offset) {
		int match[2] = {0, 0};
		for (int i = 0; i < size; ++i) {
			if (i + offset < 0) continue;
			if (i + offset >= size) continue;
			// for (int j = 0; j < 32; ++j) {
			// 	if (energy[j][i] && time[j][i+offset]) ++match[j/16];
			// }
			if (energy[0][i] && time[0][i+offset]) ++match[0];
			if (energy[16][i] && time[16][i+offset]) ++match[1];
		}
		gm[0].AddPoint(offset, match[0]);
		gm[1].AddPoint(offset, match[1]);
	}

	// find max match in GDC
	int points = gm[0].GetN();
	double *offsets = gm[0].GetX();
	double *match = gm[0].GetY();
	// max offset
	double max_offset = offsets[0];
	// max match number
	double max_match = match[0];
	for (int i = 1; i < points; ++i) {
		if (match[i] > max_match) {
			max_match = match[i];
			max_offset = offsets[i];
		}
	}
	std::cout << "GDC 0: " << max_offset << "\n";

	gm[0].Write("gm0");
	gm[1].Write("gm1");
	opf.Close();

	return 0;
}