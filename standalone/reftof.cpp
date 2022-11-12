#include <iostream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>

#include "include/defs.h"

using namespace ribll;

const unsigned int look_window = 10'000;

std::vector<double> tof_time_list;

int main(int argc, char **argv) {
	if (argc != 3) {
		std::cout << "Usage: " << argv[0] << " run detector\n";
		return -1;
	}

	int run = atoi(argv[1]);
	std::string detector = std::string(argv[2]);

	TString tof_file_name;
	tof_file_name.Form("%s%stof-map-%04d.root", kGenerateDataPath, kMappingDir, run);
	TFile *tof_file = new TFile(tof_file_name, "read");
	TTree *tof_tree = (TTree*)tof_file->Get("tree");
	if (!tof_tree) {
		std::cerr << "Error: get tree from " << tof_file_name << " failed.\n";
		return -1;
	}

	double tof2_time;
	unsigned short tof_index;
	tof_tree->SetBranchAddress("time", &tof2_time);
	tof_tree->SetBranchAddress("index", &tof_index);

	// show process
	printf("Reading tof events   0%%");
	fflush(stdout);
	// 1/100 of entry
	long long entry100 = tof_tree->GetEntries() / 100;
	for (long long entry = 0; entry < tof_tree->GetEntries(); ++entry) {
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		tof_tree->GetEntry(entry);
		if (tof_index != 1) continue;
		tof_time_list.push_back(tof2_time);
	}
	printf("\b\b\b\b100%%\n");
	tof_file->Close();

	// sort tof events
	std::sort(tof_time_list.begin(), tof_time_list.end());
	for (size_t i = 0; i < 10; ++i) {
		std::cout << tof_time_list[i] << "\n";
	}


	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-single-side-%04d.root",
		kGenerateDataPath, kSingleSideDir, detector.c_str(), run
	);
	TFile *ipf = new TFile(input_file_name, "read");
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: get tree from " << input_file_name << " failed.\n";
		return -1;
	}

	double detector_time[8];
	ipt->SetBranchAddress("time", detector_time);
	

	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-reftof-%04d.root",
		kGenerateDataPath, kSingleSideDir, detector.c_str(), run
	);
	TFile *opf = new TFile(output_file_name, "recreate");
	TH1F *hist_look_window = new TH1F("ht", "look window", 1000, -look_window, look_window);
	TTree *opt = new TTree("tree", "reference tof2");

	opt->Branch("time", &tof2_time, "t/D");

	std::cout << "ipt entries " << ipt->GetEntries() << std::endl;

	// show process
	printf("Reading detector events   0%%");
	fflush(stdout);
	// 1/100 of entry
	entry100 = ipt->GetEntries() / 100;
	for (long long entry = 0; entry < ipt->GetEntries(); ++entry) {
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);

		bool found = false;
		for (
			auto iter = std::lower_bound(tof_time_list.begin(), tof_time_list.end(), detector_time[0] - look_window);
			iter != tof_time_list.end();
			++iter
		) {
			if (*iter > detector_time[0] + look_window) break;
			hist_look_window->Fill(*iter - detector_time[0]);
			if (*iter >= detector_time[0]) continue;
			if (*iter <= detector_time[0] - 500) continue;
			if (!found) {	
				tof2_time = *iter;
				found = true;
			}
		}

		if (!found) tof2_time = -1.0;
		opt->Fill();
	}
	printf("\b\b\b\b100%%\n");
	opt->Write();
	hist_look_window->Write();
	std::cout << "opt entries " << opt->GetEntries() << std::endl;
	opf->Close();

	ipf->Close();


	return 0;
}