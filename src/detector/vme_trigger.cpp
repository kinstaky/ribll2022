#include "include/detector/vme_trigger.h"

#include "include/event/trigger_event.h"
#include "include/event/ppac_event.h"
#include "include/event/dssd_event.h"
#include "include/event/tof_event.h"

namespace ribll {

VmeTrigger::VmeTrigger(unsigned int run, const std::string &tag)
: Detector(run, "vt", tag) {
}


/// @brief convert map event to fundamental event
/// @param[in] trigger_time trigger time to match
/// @param[in] match_map map_events order by trigger time
/// @param[out] fundamental_event converted fundamental event
/// @param[inout] statistics information about statistics
/// @returns 0 if success, -1 otherwise
///
int FillEvent(
	double trigger_time,
	const std::multimap<double, TriggerEvent> &match_map,
	TriggerEvent &fundamental_event,
	std::vector<MatchTriggerStatistics> &statistics
) {
	fundamental_event.time = -1e5;
	fundamental_event.cfd_flag = 0;
	fundamental_event.timestamp = -1;

	size_t match_count = match_map.count(trigger_time);
	if (match_count == 1) {
		const auto &event = match_map.equal_range(trigger_time).first->second;
		fundamental_event.time = event.time - trigger_time;
		fundamental_event.timestamp = event.timestamp;
		fundamental_event.cfd_flag = event.cfd_flag;
		++statistics[0].match_events;
		++statistics[0].used_events;

		return 0;
	} else if (match_count > 1) {
		++statistics[0].oversize_events;
	}
	return -1;
}


template<typename FundamentalEvent>
int MatchVmeDetector(
	const std::string &trigger_tag,
	unsigned int run,
	const std::string &name
) {
	// setup reference input
	// name of reference vme trigger funamental file
	TString vme_trigger_file_name;
	vme_trigger_file_name.Form(
		"%s%svt-fundamental-%svt-%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		trigger_tag.empty() ? "" : (trigger_tag + "-").c_str(),
		run
	);
	// pointer to reference vme trigger fundamental file
	TFile *trigger_file = new TFile(vme_trigger_file_name, "read");
	// pointer to reference vme trigger fundamental tree
	TTree *trigger_tree = (TTree*)trigger_file->Get("tree");
	if (!trigger_tree) {
		std::cerr << "Error: Get tree from "
			<< vme_trigger_file_name << " failed.\n";
		return -1;
	}
	// name of XIA trigger file
	TString xia_trigger_file_name;
	xia_trigger_file_name.Form(
		"%s%sxt-map-%svt-%04u.root",
		kGenerateDataPath,
		kMappingDir,
		trigger_tag.empty() ? "" : (trigger_tag + "-").c_str(),
		run
	);
	// add friend
	trigger_tree->AddFriend("xt=tree", xia_trigger_file_name);
	// reference vme trigger time
	long long vme_trigger_time;
	// reference xia trigger time
	double xia_trigger_time;
	// setup input branch
	trigger_tree->SetBranchAddress("timestamp", &vme_trigger_time);
	trigger_tree->SetBranchAddress("xt.time", &xia_trigger_time);

	// setup detector map event input
	// name of input detector map event file
	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-map-%04u.root",
		kGenerateDataPath, kMappingDir, name.c_str(), run
	);
	// pointer to input detector map event file
	TFile *ipf = new TFile(input_file_name, "read");
	// pointer to input detector map event file
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}

	// setup detector fundamental event output
	// name of output detector fundamental event file
	TString out_detector_file_name;
	out_detector_file_name.Form(
		"%s%s%s-fundamental-%s%s-%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		name.c_str(),
		trigger_tag.empty() ? "" : (trigger_tag + "-").c_str(),
		name.c_str(),
		run
	);
	// pointer to output detector fundamental event file
	TFile *out_detector_file = new TFile(out_detector_file_name, "recreate");
	// pointer to output detector fundamental event tree
	TTree *out_detector_tree = new TTree("tree", "fundamental tree");

	// input and output event
	FundamentalEvent event;
	// input align detector align trigger time
	long long align_time;
	// setup input and output detector tree branches
	event.SetupInput(ipt);
	ipt->SetBranchAddress("timestamp", &align_time);
	event.SetupOutput(out_detector_tree);

	// setup output extracted XIA trigger file and tree
 	// name of output extracted XIA trigger file
	TString out_trigger_file_name;
	out_trigger_file_name.Form(
		"%s%sxt-map-%s%s-%04u.root",
		kGenerateDataPath,
		kMappingDir,
		trigger_tag.empty() ? "" : (trigger_tag + "-").c_str(),
		name.c_str(),
		run
	);
	// output extracted XIA trigger file
	TFile *out_trigger_file = new TFile(out_trigger_file_name, "recreate");
	// output extracted XIA trigger tree
	TTree *out_trigger_tree = new TTree("tree", "extracted trigger");
	// setup branch
	out_trigger_tree->Branch("time", &xia_trigger_time, "t/D");

	// statistics
	MatchTriggerStatistics statistics(
		run,
		name,
		trigger_tag,
		name,
		trigger_tree->GetEntries(),
		ipt->GetEntries()
	);

	// total entry of reference trigger
	long long entries = trigger_tree->GetEntries();
	// 1/100 of total entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show begin
	printf(
		"Matching %s events in VME (%s)   0%%",
		name.c_str(),
		trigger_tag.c_str()
	);
	fflush(stdout);
	// input map event entry
	long long input_entry = 0;
	// total entry of input tree
	long long total_input_entries = ipt->GetEntries();
	// loop to match reference trigger
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3llu%%", entry / entry100);
			fflush(stdout);
		}

		// get reference trigger time
		trigger_tree->GetEntry(entry);

		// reach the end of detector input file, finish
		if (input_entry >= total_input_entries) break;

		ipt->GetEntry(input_entry);
		while (
			align_time < vme_trigger_time
			&& input_entry + 1 < total_input_entries
		) {
			// If align time is smaller than trigger time, this VME trigger
			// doesn't match any XIA trigger. So continue and go to check
			// next event.
			++input_entry;
			ipt->GetEntry(input_entry);
		}
		// Now align time is larger than or equal to the VME trigger time.
		if (align_time == vme_trigger_time) {
			// Align time is equal to trigger time, which means this event
			// match the XIA trigger.
			// jump to next entry
			++input_entry;
			++statistics.match_events;
			++statistics.used_events;
			out_detector_tree->Fill();
			out_trigger_tree->Fill();
		}
		// align time may be larger than VME trigger time,
		// do nothing and continue
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// write tree and close file
	out_detector_file->cd();
	out_detector_tree->Write();
	out_detector_file->Close();

	out_trigger_file->cd();
	out_trigger_tree->Write();
	out_trigger_file->Close();

	// close reference file
	trigger_file->Close();
	// close input file
	ipf->Close();

	statistics.Write();
	statistics.Print();

	return 0;
}


int VmeTrigger::ExtractTrigger(
	double window_left,
	double window_right
) {
	if (Detector::ExtractTrigger<TriggerEvent, TriggerEvent>(
		{"vt"},
		window_left,
		window_right,
		FillEvent
	)) return -1;

	// // match detectors event recorded by VME
	// // ToF in VME
	// if (MatchVmeDetector<TofFundamentalEvent>(trigger_tag, run_, "vtof")) {
	// 	std::cerr << "Error: Match VME event in vtof failed.\n";
	// 	return -1;
	// }
	// // PPAC in VME
	// if (MatchVmeDetector<PpacFundamentalEvent>(trigger_tag, run_, "vppac")) {
	// 	std::cerr << "Error: Match VME event in vppac failed.\n";
	// 	return -1;
	// }
	// // TAF0 and TAF1 in VME
	// for (int i = 0; i < 2; ++i) {
	// 	std::string detector = "taf" + std::to_string(i);
	// 	if (MatchVmeDetector<DssdFundamentalEvent>(
	// 		trigger_tag, run_, detector
	// 	)) {
	// 		std::cerr << "Error: Match VME event in "
	// 			<< detector << " failed.\n";
	// 		return -1;
	// 	}
	// }
	// // TAB0 ~ TAB5 in VME
	// for (int i = 0; i < 6; ++i) {
	// 	std::string detector = "tab" + std::to_string(i);
	// 	if (MatchVmeDetector<DssdFundamentalEvent>(
	// 		trigger_tag, run_, detector
	// 	)) {
	// 		std::cerr << "Error: Match VME event in "
	// 			<< detector << " failed.\n";
	// 		return -1;
	// 	}
	// }

	return 0;
}


}