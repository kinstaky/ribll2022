#include "include/detector/ppac.h"

#include <iostream>
#include <map>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

#include "include/defs.h"

namespace ribll {

const long long search_window = 1000;
const long long look_window = 10000; 

PPAC::PPAC(int run)
:run_(run) {
}

// int PPAC::ReadTrigger() {
// 	events_.clear();
// 	// setup input file and tree
// 	TString file_name;
// 	file_name.Form("%s%sxt-map-%04d.root", kGenerateDataPath, kMappingDir, run_);
// 	TFile *ipf = new TFile(file_name, "read");
// 	if (!ipf) {
// 		std::cerr << "Error: open file " << file_name << " failed.\n";
// 		return -1;
// 	}
// 	TTree *ipt = (TTree*)ipf->Get("tree");
// 	if (!ipt) {
// 		std::cerr << "Error: get tree from " << file_name << " failed.\n";
// 		ipf->Close();
// 		return -1;
// 	}

// 	// input data and branch
// 	Long64_t timestamp;


// 	ipf->Close();
// }


int PPAC::ReadEvents() {
	events_.clear();
	// setup input file and tree
	TString file_name;
	file_name.Form("%s%sxppac-map-%04d.root", kGenerateDataPath, kMappingDir, run_);
	TFile *ipf = new TFile(file_name, "read");
	if (!ipf) {
		std::cerr << "Error: open file " << file_name << " failed.\n";
		return -1;
	}
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: get tree from " << file_name << " failed.\n";
		ipf->Close();
		return -1;
	}

	// input data
	Long64_t timestamp;
	// setup input branch
	ipt->SetBranchAddress("timestamp", &timestamp);
	ipt->SetBranchAddress("index", &event_.index);
	ipt->SetBranchAddress("side", &event_.side);
	ipt->SetBranchAddress("time", &event_.time);
	event_.used = false;
	
	// show process
	printf("Reading ppac events   0%%");
	fflush(stdout);
	Long64_t entry100 = ipt->GetEntries() / 100;
	for (Long64_t entry = 0; entry < ipt->GetEntries(); ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);
		events_.insert(std::make_pair(timestamp, event_));
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	return 0;
}

int PPAC::Correlate() {
	// read ppac events
	if (ReadEvents()) {
		std::cerr << "Error: read ppac events failed.\n";
		return -1;
	}

	// setup input file
	TString input_file_name;
	input_file_name.Form("%s%sxt-map-%04d.root", kGenerateDataPath, kMappingDir, run_);
	TFile *ipf = new TFile(input_file_name, "read");
	if (!ipf) {
		std::cerr << "Error: open file " << input_file_name << " failed.\n";
		return -1;
	}
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: get tree from " << input_file_name << " failed.\n";
		ipf->Close();
		return -1;
	}
	// input data and branch
	Long64_t timestamp;
	// setup input trigger branch
	ipt->SetBranchAddress("timestamp", &timestamp);            


	// setup output file
	TString output_file_name;
	output_file_name.Form("%s%sppac-corr-%04d.root", kGenerateDataPath, kCorrelationDir, run_); 
	TFile *opf = new TFile(output_file_name, "recreate");
	if (!opf) {
		std::cerr << "Error: create file " << output_file_name << " failed.\n";
		return -1;
	}
	TH1F *hist_look_window = new TH1F(
		"ht", "time window", 1000, -look_window, look_window
	);
	TTree *opt = new TTree("tree", "correlation tree of ppac");
	// output data
	long long ppac_timestamp;
	int ppac_flag;
	unsigned short ppac_hit;
	unsigned short ppac_xhit;
	unsigned short ppac_yhit;
	double ppac_x[ppac_num];
	double ppac_y[ppac_num];
	// setup output branches
	opt->Branch("timestamp", &ppac_timestamp, "ts/L");
	opt->Branch("flag", &ppac_flag, "flag/I");
	opt->Branch("hit", &ppac_hit, "hit/s");
	opt->Branch("xhit", &ppac_xhit, "xhit/s");
	opt->Branch("yhit", &ppac_yhit, "yhit/s");
	opt->Branch("x", ppac_x, TString::Format("x[%llu]/D", ppac_num).Data());
	opt->Branch("y", ppac_y, TString::Format("y[%llu]/D", ppac_num).Data());


	// show process
	printf("Correlating ppac events   0%%");
	fflush(stdout);
	Long64_t entry100 = ipt->GetEntries() / 100;
	long long correlation_num = 0;
	for (Long64_t entry = 0; entry < ipt->GetEntries(); ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get trigger timestamp
		ipt->GetEntry(entry);

		std::vector<Event> filling_events;
		
		// search events
		ppac_flag = 0;
		ppac_hit = ppac_xhit = ppac_yhit = 0;
		for (size_t i = 0; i < ppac_num; ++i) {
			ppac_x[i] = 0.0;
			ppac_y[i] = 0.0;
		}
		for (
			auto iter = events_.lower_bound(timestamp - look_window);
			iter != events_.upper_bound(timestamp + look_window);
			++iter
		) {
			hist_look_window->Fill(iter->first - timestamp);
			
			if (iter->first <= timestamp - search_window) continue;
			if (iter->first >= timestamp + search_window) continue;
			if (iter->second.used) continue;
			int flag_bit = 1 << (iter->second.index * 5 + iter->second.side);
			if ((flag_bit & ppac_flag) != 0) continue;

			iter->second.used = true;
			if (!ppac_flag) {
				// first event
				ppac_timestamp = iter->first;
			}
			ppac_flag |= flag_bit;
			++ppac_hit;
			filling_events.push_back(iter->second);
		}

		// sort events
		std::sort(filling_events.begin(), filling_events.end(), [](const Event &x, const Event &y) {
			return (x.index < y.index) || (x.index == y.index && x.side < y.side);
		});

		// fill to tree
		if (filling_events.empty()) continue;
		for (size_t i = 0; i < filling_events.size()-1; ++i) {
			if (filling_events[i].side == 0 && filling_events[i+1].side == 1) {
				int bit_flag
					= 3 << (filling_events[i].index * 5 + filling_events[i].side);
				if ((bit_flag & ppac_flag) == bit_flag) {
					++ppac_xhit;
					ppac_x[filling_events[i].index]
						= filling_events[i].time - filling_events[i+1].time;
				}
			} else if (filling_events[i].side == 2 && filling_events[i+1].side == 3) {
				int bit_flag
					= 3 << (filling_events[i].index * 5 + filling_events[i].side);
				if ((bit_flag & ppac_flag) == bit_flag) {
					++ppac_yhit;
					ppac_y[filling_events[i].index]
						= filling_events[i].time - filling_events[i+1].time;
				}
			}
		}
		if (ppac_xhit && ppac_yhit) {
			opt->Fill();
			correlation_num++;
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	hist_look_window->Write();
	opt->Write();
	opf->Close();

	// show correlation rate
	std::cout << "ppac correlation rate " << correlation_num << " / " << ipt->GetEntries()
		<< " " << double(correlation_num) / double(ipt->GetEntries()) << "\n";


	ipf->Close();

	
	return 0;
}


}