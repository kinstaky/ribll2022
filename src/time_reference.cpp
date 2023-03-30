#include "include/time_reference.h"

#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

namespace ribll {

int TimeReference(unsigned int run, const std::string &tag) {
	// input xppac file
	TString ppac_file_name;
	ppac_file_name.Form(
		"%s%sxppac-fundamental-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// ppac file
	TFile *ppac_file = new TFile(ppac_file_name, "read");
	// ppac tree
	TTree *ppac_tree = (TTree*)ppac_file->Get("tree");
	if (!ppac_tree) {
		std::cerr << "Error: Get tree from "
			<< ppac_file_name << " failed.\n";
		return -1;
	}
	// input xppac event
	PpacFundamentalEvent ppac_event;
	// setup input branches
	ppac_event.SetupInput(ppac_tree);

	// output reference file name
	TString ref_file_name;
	ref_file_name.Form(
		"%s%sreftime-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// output reference file
	TFile *ref_file = new TFile(ref_file_name, "recreate");
	// output reference tree
	TTree *ref_tree = new TTree("tree", "reference time");
	// output reference time
	double ref_time;
	// setup output tree branch
	ref_tree->Branch("time", &ref_time, "t/D");

	// total number of entries
	long long entries = ppac_tree->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Calculating reference time   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// read ppac event
		ppac_tree->GetEntry(entry);
		// initialize
		ref_time = -1e5;
		if (
			(ppac_event.flag & 0x4200) == 0x4200
			&& ppac_event.anode[1] - ppac_event.anode[2] > -16
			&& ppac_event.anode[1] - ppac_event.anode[2] < -13
		) {
			// a1 and a2
			ref_time = ppac_event.anode[2];
		} else if (
			(ppac_event.flag & 0x4010) == 0x4010
			&& ppac_event.anode[0] - ppac_event.anode[2] > -18
			&& ppac_event.anode[0] - ppac_event.anode[2] < -13
		) {
			// a0 and a2
			ref_time = ppac_event.anode[2];
		} else if (
			(ppac_event.flag & 0x210) == 0x210
			&& ppac_event.anode[0] - ppac_event.anode[1] > -2
			&& ppac_event.anode[0] - ppac_event.anode[1] < 1.5		
		) {
			// a0 and a1
			ref_time = ppac_event.anode[1] + 14.2;
		}
		ref_tree->Fill();
	}
	printf("\b\b\b\b100%%\n");

	// save tree and close files
	ref_tree->Write();
	ref_file->Close();
	ppac_file->Close();

	return 0;
}

}		// time reference