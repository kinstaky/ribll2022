#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TRandom3.h>

#include "include/event/generate_event.h"
#include "include/event/tof_event.h"

using namespace ribll;

constexpr double vtof_time[2] = {-522.0, -462.0};

int main(int argc, char **argv) {
	int run = 0;
	if (argc > 1) {
		run = atoi(argv[1]);
	}
	if (run < 0 || run > 4) {
		std::cout << "Usage: " << argv[0] << "[run]\n"
			<< "  run        run number, default is 0\n";
		return -1;
	}

	// input generate file name
	TString input_file_name = TString::Format(
		"%s%sgenerate-%04d.root", kGenerateDataPath, kSimulateDir, run
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
	// input generate event
	GenerateEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// output file name
	TString output_file_name = TString::Format(
		"%s%svtof-fundamental-sim-ta-%04d.root",
		kGenerateDataPath,
		kFundamentalDir,
		run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "simulated VME ToF fundamental events");
	// output event
	TofFundamentalEvent vtof;
	//setup output branches
	vtof.SetupOutput(&opt);

	// initialize
	vtof.cfd_flag = 0;
	// random number generator
	TRandom3 generator(0);

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100;
	// show start
	printf("Simulating detect VME ToF   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get data
		ipt->GetEntry(entry);

		vtof.time[0] = generator.Gaus(vtof_time[0], 0.15);
		vtof.time[1] = generator.Gaus(vtof_time[1], 0.15);

		// fill tree
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}