#include <iostream>
#include <map>

#include <TH1F.h>

#include "include/detector/tof.h"
#include "include/defs.h"
#include "include/event/tof_event.h"


namespace ribll {

const double look_window = 10000;

Tof::Tof(int run)
: run_(run) {
}

int Tof::ReadTriggerTimes(std::vector<double> &trigger_times) {
	trigger_times.clear();
	// trigger file name
	TString trigger_file_name;
	trigger_file_name.Form(
		"%s%sxt-map-%04d.root",
		kGenerateDataPath, kMappingDir, run_
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
	printf("Reading trigger events   0%%");
	fflush(stdout);
	long long entries = trigger_tree->GetEntries();
	// 1/100 of entry, for showing process
	long long entry100 = entries / 100;
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


int Tof::MatchTrigger(double window_left, double window_right) {
	std::vector<double> trigger_times;
	if (ReadTriggerTimes(trigger_times)) {
		return -1;
	}

	// setup input file and tree
	// name of input file
	TString input_file_name;
	input_file_name.Form(
		"%s%stof-map-%04d.root",
		kGenerateDataPath, kMappingDir, run_
	);
	// pointer to input file
	TFile *ipf = new TFile(input_file_name, "read");
	// pointer to input tree
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// setup branches
	TofMapEvent map_event;
	map_event.SetupInput(ipt);


	// setup output file and tree
	// name of output file
	TString output_file_name;
	output_file_name.Form(
		"%s%stof-fundamental-%04d.root",
		kGenerateDataPath, kFundamentalDir, run_
	);
	// output file
	TFile *opf = new TFile(output_file_name, "recreate");
	// histogram of look window
	TH1F *hist_look_window[2];
	hist_look_window[0] = new TH1F("ht1", "look window", 1000, -look_window, look_window);
	hist_look_window[1] = new TH1F("ht2", "look window", 1000, -look_window, look_window);
	// output tree
	TTree *opt = new TTree("tree", "reference tree");
	// reference time found from detector
	TofFundamentalEvent fundamental_event;
	fundamental_event.SetupOutput(opt);

	std::multimap<double, TofMapEvent> match_map;

	// show begin
	printf("Matching detector events   0%%");
	fflush(stdout);
	// 1/100 of entry
	long long entries = ipt->GetEntries();
	long long entry100 = entries / 100;
	// loop for matching trigger time
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get entry
		ipt->GetEntry(entry);

		// 0 for not found
		double match_time = 0.0;
		for (
			auto iter = std::lower_bound(
				trigger_times.begin(),
				trigger_times.end(),
				map_event.time - look_window
			);
			iter != trigger_times.end();
			++iter
		) {
			// stop if out of look window
			if (*iter > map_event.time + look_window) break;
			// add to histogram for later check
			hist_look_window[map_event.index]->Fill(*iter - map_event.time);
			// jump if out of search window
			if (*iter <= map_event.time + window_left) continue;
			if (*iter >= map_event.time + window_right) continue;

			if (
				fabs(*iter-map_event.time) < fabs(match_time-map_event.time)
			) {
				// found the closer referece trigger time
				match_time = *iter;
			}
		}

		if (match_time != 0.0) {
			// increment if found reference events for stastics
			match_map.insert(std::make_pair(match_time, map_event));
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// write the histogram
	hist_look_window[0]->Write();
	hist_look_window[1]->Write();
	// close input file
	ipf->Close();


	// for statistics
	long long match_events = 0;

	// show begin
	printf("Writing fundamental events   0%%");
	fflush(stdout);
	entries = trigger_times.size();
	entry100 = entries / 100;
	// loop for matching trigger time
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// initialize output event
		fundamental_event.time[0] = -1.0;
		fundamental_event.time[1] = -1.0;

		// check match events number
		size_t match_count = match_map.count(trigger_times[entry]);
		if (match_count == 1 || match_count == 2) {
			auto range = match_map.equal_range(trigger_times[entry]);
			for (auto iter = range.first; iter != range.second; ++iter) {
				fundamental_event.time[iter->second.index] = iter->second.time;
			}
			match_events++;
		}
		opt->Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// write the tree
	opt->Write();
	// close output file
	opf->Close();

	// show statistics
	std::cout << "events match rate " << match_events << " / " << entries
		<< "  " << double(match_events) / double(entries) << "\n";

	return 0;
}

}