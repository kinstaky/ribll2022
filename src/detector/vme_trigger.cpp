#include "include/detector/vme_trigger.h"

#include "include/event/trigger_event.h"
#include "include/event/ppac_event.h"
#include "include/event/dssd_event.h"

namespace ribll {

VmeTrigger::VmeTrigger(unsigned int run)
: Detector(run, "vt") {
}


/// @brief convert map event to fundamental event
/// @param[in] trigger_time trigger time to match
/// @param[in] match_map map_events order by trigger time
/// @param[out] fundamental_event converted fundamental event
/// @param[inout] statistics information about statistics
///
void FillEvent(
	double trigger_time,
	const std::multimap<double, TriggerEvent> &match_map,
	TriggerEvent &fundamental_event,
	MatchTriggerStatistics &statistics
) {
	fundamental_event.time = -1e5;
	fundamental_event.cfd_flag = 0;
	fundamental_event.timestamp = -1;

	size_t match_count = match_map.count(trigger_time);
	if (match_count == 1) {
		auto range = match_map.equal_range(trigger_time);
		for (auto iter = range.first; iter != range.second; ++iter) {
			fundamental_event.time = iter->second.time - trigger_time;
			fundamental_event.timestamp = iter->second.timestamp;
			fundamental_event.cfd_flag = iter->second.cfd_flag;
		}
		++statistics.match_events;
		++statistics.used_events;
	} else if (match_count > 1) {
		++statistics.oversize_events;
	}
}


template<typename FundamentalEvent>
int MatchVmeDetector(unsigned int run, const std::string &name) {
	// setup reference input
	// name of reference vme trigger funamental file 
	TString vme_trigger_file_name;
	vme_trigger_file_name.Form(
		"%s%svt-fundamental-%04u.root",
		kGenerateDataPath, kFundamentalDir, run
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
	// reference vme trigger time
	long long trigger_time;
	// setup input branch
	trigger_tree->SetBranchAddress("timestamp", &trigger_time);

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
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-fundamental-%04u.root",
		kGenerateDataPath, kFundamentalDir, name.c_str(), run
	);
	// pointer to output detector fundamental event file
	TFile *opf = new TFile(output_file_name, "recreate");
	// pointer to output detector fundamental event tree
	TTree *opt = new TTree("tree", "fundamental tree");

	// input and output event
	FundamentalEvent event;
	// input align detector align trigger time
	long long align_time;
	// setup input and output detector tree branches
	event.SetupInput(ipt);
	ipt->SetBranchAddress("time", &align_time);
	event.SetupOutput(opt);

	MatchTriggerStatistics statistics(
		run,
		name,
		trigger_tree->GetEntries(),
		ipt->GetEntries()
	);

	// total entry of reference trigger
	long long entries = trigger_tree->GetEntries();
	// 1/100 of total entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show begin
	printf("Matching %s events in VME   0%%", name.c_str());
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

		if (trigger_time < 0) {
			// reference trigger not found, impossible to match VME events
			event.Nullify();
		} else if (input_entry >= total_input_entries) {
			// reach the end of input file
			event.Nullify();
		} else {
			// reference trigger found, check input event
			ipt->GetEntry(input_entry);
			while (
				align_time < trigger_time
				&& input_entry + 1 < total_input_entries
			) {
				// If align time is smaller than trigger time, this VME trigger
				// doesn't match any XIA trigger. So continue and go to check
				// next event.
				++input_entry;
				ipt->GetEntry(input_entry);
			}
			// Now align time is larger than or equal to the trigger time.
			if (align_time == trigger_time) {
				// Align time is equal to trigger time, which means this event
				// match the XIA trigger.
				// jump to next entry
				++input_entry;
				++statistics.match_events;
				++statistics.used_events;
			} else {
				// Align time is larger than trigger time, which mean no
				// match events is found.
				event.Nullify();
			}
		}

		opt->Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// write tree
	opt->Write();

	// close reference file
	trigger_file->Close();
	// close input file
	ipf->Close();
	// close output file
	opf->Close();

	statistics.Write();
	statistics.Print();

	return 0;
}


int VmeTrigger::MatchTrigger(double window_left, double window_right) {
	if (Detector::MatchTrigger<TriggerEvent, TriggerEvent>(
		window_left,
		window_right,
		FillEvent
	)) {
		return -1;
	}

	// match detectors event recorded by VME
	// PPAC in VME
	if (MatchVmeDetector<PpacFundamentalEvent>(run_, "vppac")) {
		std::cerr << "Error: Match VME event in vppac failed.\n";
		return -1;
	}
	// TAF0 and TAF1 in VME
	for (int i = 0; i < 2; ++i) {
		std::string detector = "taf" + std::to_string(i);
		if (MatchVmeDetector<DssdFundamentalEvent>(run_, detector)) {
			std::cerr << "Error: Match VME event in "
				<< detector << " failed.\n";
			return -1;
		}
	}
	// TAB0 ~ TAB5 in VME
	for (int i = 0; i < 6; ++i) {
		std::string detector = "tab" + std::to_string(i);
		if (MatchVmeDetector<DssdFundamentalEvent>(run_, detector)) {
			std::cerr << "Error: Match VME event in "
				<< detector << " failed.\n";
			return -1;
		}
	}

	return 0;
}


}