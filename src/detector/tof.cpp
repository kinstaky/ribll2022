#include "include/detector/tof.h"

#include <iostream>
#include <map>

#include <TH1F.h>
#include <TString.h>
#include <TF1.h>

#include "include/defs.h"
#include "include/event/tof_event.h"
#include "include/statistics/beam_identify_statistics.h"


namespace ribll {


Tof::Tof(unsigned int run)
: Detector(run, "tof") {
}


// /// @brief convert map event to fundamental event
// /// @param[in] trigger_time trigger time to match
// /// @param[in] match_map map_events order by trigger time
// /// @param[out] fundamental_event converted fundamental event
// /// @param[out] statistics information about statistics
// /// @returns index >= 0 if matched, -1 for failed
// ///
// int FillEvent(
// 	double trigger_time,
// 	const std::multimap<double, TofMapEvent> &match_map,
// 	TofFundamentalEvent &fundamental_event,
// 	std::vector<MatchTriggerStatistics> &statistics
// ) {
// 	fundamental_event.time[0] = -1e5;
// 	fundamental_event.time[1] = -1e5;
// 	fundamental_event.cfd_flag = 0;

// 	// check match events number
// 	size_t match_count = match_map.count(trigger_time);

// 	// jump if match not found
// 	if (match_count == 0) return -1;
// 	// record oversize events and jump
// 	if (match_count > 2) {
// 		for (MatchTriggerStatistics &sta : statistics) {
// 			++sta.oversize_events;
// 		}
// 		return -1;
// 	}

// 	// index is 0 represents golden if match count is 2,
// 	// index is 1 represents silver if match count is 1
// 	size_t index = 2 - match_count;
// 	auto range = match_map.equal_range(trigger_time);
// 	for (auto iter = range.first; iter != range.second; ++iter) {
// 		fundamental_event.time[iter->second.index] =
// 			iter->second.time - trigger_time;
// 		fundamental_event.cfd_flag |=
// 			iter->second.cfd_flag ? (1 << iter->second.index) : 0;
// 	}
// 	++statistics[index].match_events;
// 	statistics[index].used_events += match_count;

// 	return index;
// }


// int Tof::ExtractTrigger(
// 	const std::string &trigger_tag,
// 	double window_left,
// 	double window_right
// ) {
// 	return Detector::ExtractTrigger<TofMapEvent, TofFundamentalEvent>(
// 		trigger_tag,
// 		{"golden", "silver"},
// 		window_left,
// 		window_right,
// 		FillEvent
// 	);
// }


/// @brief convert map event to fundamental event
/// @param[in] trigger_time trigger time to match
/// @param[in] match_map map_events order by trigger time
/// @param[out] fundamental_event converted fundamental event
/// @param[out] statistics information about statistics
///
void FillEvent(
	double trigger_time,
	const std::multimap<double, TofMapEvent> &match_map,
	TofFundamentalEvent &fundamental_event,
	MatchTriggerStatistics &statistics
) {
	fundamental_event.time[0] = -1e5;
	fundamental_event.time[1] = -1e5;
	fundamental_event.cfd_flag = 0;

	// check match events number
	size_t match_count = match_map.count(trigger_time);

	if (match_count == 1 || match_count == 2) {
		auto range = match_map.equal_range(trigger_time);
		for (auto iter = range.first; iter != range.second; ++iter) {
			fundamental_event.time[iter->second.index] =
				iter->second.time - trigger_time;
			fundamental_event.cfd_flag |=
				iter->second.cfd_flag ? (1 << iter->second.index) : 0;
		}
		++statistics.match_events;
		statistics.used_events += match_count;
	} else if (match_count > 2) {
		++statistics.oversize_events;
	}
}


int Tof::MatchTrigger(
	const std::string &trigger_tag,
	double window_left,
	double window_right
) {
	if (name_ == "tof") {
		return Detector::MatchTrigger<TofMapEvent, TofFundamentalEvent>(
			trigger_tag,
			window_left,
			window_right,
			FillEvent
		);
	} else {
		// vtof
		return Detector::VmeMatchTrigger<TofFundamentalEvent>(
			trigger_tag
		);
	}
}




int Tof::BeamIdentify() {
	// setup input file
	// input file name
	TString input_file_name;
	input_file_name.Form(
		"%s%stof-fundamental-%04u.root",
		kGenerateDataPath, kFundamentalDir, run_
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
	// input event
	TofFundamentalEvent fundamental_event;
	// setup input branches
	fundamental_event.SetupInput(ipt);

	// setup output file
	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%sbeam-type-%04u.root",
		kGenerateDataPath, kBeamDir, run_
	);
	// output file
	TFile *opf = new TFile(output_file_name, "recreate");
	// output ToF histogram and fitting
	TH1F *htof = new TH1F("htof", "ToF(T2-T1) of beam", 1000, 0, 100);
	// output tree
	TTree *opt = new TTree("tree", "tree of beam type");
	// output beam type
	short beam_type;
	// output branches
	opt->Branch("beam_type", &beam_type, "bt/S");

	// first loop to fill histogram
	// total number of input tree
	long long entries = ipt->GetEntries();
	for (long long entry = 0; entry < entries; ++entry) {
		ipt->GetEntry(entry);

		// ignoe invalid tof value
		if (
			fundamental_event.time[0] < -9e4
			|| fundamental_event.time[1] < -9e4
		) continue;

		htof->Fill(fundamental_event.time[1] - fundamental_event.time[0]);
	}

	// fit peaks
	TF1 *f14c = new TF1("f14c", "gaus", 50, 70);
	f14c->SetParameter(0, 1e5);
	f14c->SetParameter(1, 60);
	f14c->SetParameter(2, 2);
	f14c->SetParLimits(1, 50, 70);
	htof->Fit(f14c, "BQR+");
	// get fit result
	BeamIdentifyStatistics statistics(run_, entries);

	// constant of 14C
	statistics.const14c = f14c->GetParameter(0);
	// mean of 14C
	statistics.mean14c = f14c->GetParameter(1);
	// sigma of 14C
	statistics.sigma14c = f14c->GetParameter(2);

	// save histograms with fitting
	htof->Write();

	// second loop to record beam type
	// minimum tof value of 14C
	double min14c = statistics.mean14c - 5 * statistics.sigma14c;
	double max14c = statistics.mean14c + 5 * statistics.sigma14c;
	for (long long entry = 0; entry < entries; ++entry) {
		ipt->GetEntry(entry);

		// initialize beam type
		beam_type = 0;

		// ignoe invalid tof value
		if (
			fundamental_event.time[0] < -9e4
			|| fundamental_event.time[1] < -9e4
		) {
			beam_type = -1;
			opt->Fill();
			continue;
		}

		// tof of this particle
		double tof = fundamental_event.time[1] - fundamental_event.time[0];

		// check if it's 14C
		if (tof > min14c && tof < max14c) {
			beam_type = 1;
			++statistics.c14;
		}
		opt->Fill();
	}
	opt->Write();
	// close files
	opf->Close();
	ipf->Close();


	statistics.Write();
	statistics.Print();

	return 0;
}


}		// namespace ribll