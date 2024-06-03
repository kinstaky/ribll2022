#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TH1F.h>

#include "include/event/dssd_event.h"

using namespace ribll;

int main(int argc, char **argv) {
	if (argc != 3) {
		std::cout << "Usage: " << argv[0] << " run DSSD\n"
			<< "  run          run number\n"
			<< "  DSSD         DSSD detector name\n";
		return -1;
	}

	int run = atoi(argv[1]);
	std::string dssd(argv[2]);

	// input file name
	TString input_file_name = TString::Format(
		"%s%s%s-result-%04d.root",
		kGenerateDataPath, kNormalizeDir, dssd.c_str(), run
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
	// input event
	DssdFundamentalEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// output file name
	TString output_file_name = TString::Format(
		"%s%s%s-pixel-energy-%04d.root",
		kGenerateDataPath, kShowDir, dssd.c_str(), run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output histograms
	int strips = dssd == "t0d1" ? 32 : 16;
	int start_strip = dssd == "t0d1" ? 16 : 8;
	std::vector<TH1F> hist_front_energy;
	std::vector<TH1F> hist_back_energy;
	std::vector<TH1F> hist_energy_diff;
	for (int i = start_strip; i < start_strip+strips; ++i) {
		for (int j = start_strip; j < start_strip+strips; ++j) {
			hist_front_energy.emplace_back(
				TString::Format("hfef%db%d", i, j),
				TString::Format("front energy of pixel fs %d, bs %d", i, j),
				1000, 1000, 41000
			);
			hist_back_energy.emplace_back(
				TString::Format("hbef%db%d", i, j),
				TString::Format("back energy of pixel fs %d, bs %d", i, j),
				1000, 1000, 41000
			);
			hist_energy_diff.emplace_back(
				TString::Format("hdef%db%d", i, j),
				TString::Format("energy difference of fs %d, bs %d", i, j),
				1000, -5000, 5000
			);
		}
	}


	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling pixel energy   0%%");
	fflush(stdout);
	// loop
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get data
		ipt->GetEntry(entry);

		for (int i = 0; i < event.front_hit; ++i) {
			if (event.front_strip[i] < start_strip) continue;
			if (event.front_strip[i] >= start_strip+strips) continue;
			for (int j = 0; j < event.back_hit; ++j) {
				if (event.back_strip[j] < start_strip) continue;
				if (event.back_strip[j] >= start_strip+strips) continue;
				int index =
					(event.front_strip[i]-start_strip)*strips
					+ (event.back_strip[j]-start_strip);
				hist_front_energy[index].Fill(event.front_energy[i]);
				hist_back_energy[index].Fill(event.back_energy[j]);
				hist_energy_diff[index].Fill(
					event.front_energy[i] - event.back_energy[j]
				);
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save histograms
	for (auto &hist : hist_front_energy) hist.Write();
	for (auto &hist : hist_back_energy) hist.Write();
	for (auto &hist : hist_energy_diff) hist.Write();
	// close file
	opf.Close();
	ipf.Close();

	return 0;
}