#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

#include "include/defs.h"
#include "include/statistics.h"

using namespace ribll;

int main(int argc, char **argv) {

	// check arguments
	if (argc != 2) {
		std::cout << "Usage: " << argv[0] << " run \n"
			<< "  run                 run number\n";
		return -1;
	}
	// run number
	unsigned int run = atoi(argv[1]);

	// file name
	TString input_file_name;
	input_file_name.Form(
		"%s%sxt-map-%04u.root",
		kGenerateDataPath, kMappingDir, run
	);
	// input file
	TFile *ipf = new TFile(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// input trigger time
	double trigger_time;
	// setup branches
	ipt->SetBranchAddress("time", &trigger_time);

	TString output_file_name;
	output_file_name.Form(
		"%s%sxt-period-%04u.root",
		kGenerateDataPath, kShowDir, run
	);
	// output file
	TFile *opf = new TFile(output_file_name, "recreate");
	// output histogram
	TH1F *hp = new TH1F("hp", "time period of XIA trigger", 1000, 0, 100'000);

	// last entry's trigger time
	double last_time;
	// get the first entry and fill to last time variable
	ipt->GetEntry(0);
	last_time = trigger_time;

	// longer than one century
	double min_time = 1e20;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of total entries, for showing process
	long long entry100 = entries / 100;
	// show start
	std::cout << "Filling histogram   0%";
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
		// period
		double period = trigger_time - last_time;
		// fill period
		hp->Fill(period);
		// check if shorter than min time
		min_time = period < min_time ? period : min_time;
		// update trigger time of last entry
		last_time = trigger_time;
	}
	// show end
	std::cout << "\b\b\b\b100%\n";

	// write histogram to output file
	hp->Write();
	// close files
	opf->Close();
	ipf->Close();

	XiaTriggerPeriodStatistics statistics(run, min_time);
	statistics.Write();
	statistics.Print();

	return 0;
}