#include "include/detector/detector.h"

#include <iostream>
#include <fstream>
#include <map>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TChain.h>
#include <TGraph.h>
#include <TF1.h>

#include "include/defs.h"


namespace ribll {

Detector::Detector(
	unsigned int run,
	const std::string &name,
	const std::string &tag
)
: run_(run)
, name_(name)
, tag_(tag) {
}

int Detector::ReadTriggerTimes(std::vector<double> &trigger_times) {
	// clear data
	trigger_times.clear();
	// trigger file name
	TString trigger_file_name;
	trigger_file_name.Form(
		"%s%sxt-map-%s%04d.root",
		kGenerateDataPath,
		kMappingDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// pointer to trigger file
	TFile *trigger_file = new TFile(trigger_file_name, "read");
	// pointer to trigger tree
	TTree *trigger_tree = (TTree*)trigger_file->Get("tree");
	if (!trigger_tree) {
		std::cerr << "Error: get tree from "
			<< trigger_file_name << " failed.\n";
		return -1;
	}

	// trigger time
	double trigger_time;
	trigger_tree->SetBranchAddress("time", &trigger_time);

	// show begin
	printf("Reading %s trigger events   0%%", tag_.c_str());
	fflush(stdout);
	long long entries = trigger_tree->GetEntries();
	// 1/100 of entry, for showing process
	long long entry100 = entries / 100 + 1;
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get event
		trigger_tree->GetEntry(entry);
		trigger_times.push_back(trigger_time);
	}
	// show end
	printf("\b\b\b\b100%%\n");
	// close file
	trigger_file->Close();
	return 0;
}


int Detector::MatchTrigger(double, double) {
	// do nothing but report error
	std::cerr << "Error: MatchTrigger is not implemented yet.\n";
	return -1;
}


int Detector::ExtractTrigger(double, double) {
	// do nothing but report error
	std::cerr << "Error: ExtractTrigger is not implemented yet.\n";
	return -1;
}


}		// namespace ribll