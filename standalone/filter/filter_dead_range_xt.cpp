#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "include/event/trigger_event.h"

using namespace ribll;

int main(int argc, char **argv) {
	if (argc != 2) {
		std::cout << "Usage: " << argv[0] << " run\n"
			<< "  run              run number.\n";
		return -1;
	}
	// run number
	int run = atoi(argv[1]);

	// in microseconds
	const double dead_time = 8'640;

	// file name
	TString input_file_name;
	input_file_name.Form(
		"%s%sxt-map-%04d.root",
		kGenerateDataPath, kMappingDir, run
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// input trigger event
	TriggerEvent event;
	// setup branches
	event.SetupInput(ipt);

	TString output_file_name;
	output_file_name.Form(
		"%s%sxt-map-dead-%04u.root",
		kGenerateDataPath, kMappingDir, run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "dead time range trigger");
	// setup output branches
	event.SetupOutput(&opt);

	// last entry's trigger time
	double last_time;
	// get the first entry and fill to last time variable
	ipt->GetEntry(0);
	last_time = event.time;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of total entries, for showing process
	long long entry100 = entries / 100;
	// show start
	std::cout << "Filtering events   0%";
	fflush(stdout);
	// loop to fill histogram
	for (long long entry = 1; entry < entries; ++entry) {
		// showing process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get entry
		ipt->GetEntry(entry);
		// inverval
		double interval = event.time - last_time;
		// update trigger time of last entry
		last_time = event.time;
		if (interval < dead_time) {
			opt.Fill();
		}
	}
	// show end
	std::cout << "\b\b\b\b100%\n";

	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();

	return 0;
}