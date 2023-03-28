#ifndef __DETECTOR_H__
#define __DETECTOR_H__

#include <vector>
#include <iostream>

#include <TFile.h>
#include <TH1F.h>
#include <TString.h>
#include <TTree.h>

#include "include/defs.h"
#include "include/event/trigger_event.h"
#include "include/statistics/match_trigger_statistics.h"

namespace ribll {

class Detector {
public:

	/// @brief constructor
	/// @param[in] run run number
	/// @param[in] name detector name
	///
	Detector(unsigned int run, const std::string &name);


	/// @brief default destructor
	///
	virtual ~Detector() = default;

	//-------------------------------------------------------------------------
	//							MatchTrigger
	//-------------------------------------------------------------------------

	/// @brief match xia main trigger and build events
	/// @param[in] trigger_tag tag of trigger to chosse file
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int MatchTrigger(
		const std::string &trigger_tag,
		double window_left,
		double window_right
	);


	/// @brief template of match trigger
	/// @tparam MapEvent type of map event
	/// @tparam FundamentalEvent type of fundamental event
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @param[in] trigger_tag tag of trigger for selecting file
	/// @param[in] fill_event function to convert map event
	///		to fundamental event
	/// @returns 0 if success, -1 otherwise
	///
	template<typename MapEvent, typename FundamentalEvent>
	int MatchTrigger(
		const std::string &trigger_tag,
		double window_left,
		double window_right,
		void (*fill_event)(
			double trigger_time,
			const std::multimap<double, MapEvent> &match_map,
			FundamentalEvent &fundamental_event,
			MatchTriggerStatistics &statistics
		)
	);


	/// @brief template for VME detectors matching trigger
	/// @tparam FundamentalEvent type of fundamental event
	/// @param[in] trigger_tag tag for trigger
	/// @returns 0 if success, -1 otherwise
	/// 
	template<typename FundamentalEvent>
	int VmeMatchTrigger(const std::string &trigger_tag);


	/// @brief extract trigger with detector events
	/// @param[in] trigger_tag extract from trigger with this tag
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int ExtractTrigger(
		const std::string &trigger_tag,
		double window_left,
		double window_right
	);


	/// @brief template of extract trigger
	/// @tparam MapEvent type of map event
	/// @tparam FundamentalEvent type of fndamental event
	/// @param[in] trigger_tag extract from trigger with this tag
	/// @param[in] extract_tags tags after extraction
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @param[in] fill_event function to convert map event
	/// @returns 0 if success, -1 otherwise
	///
	template<typename MapEvent, typename FundamentalEvent>
	int ExtractTrigger(
		const std::string &trigger_tag,
		const std::vector<std::string> &extract_tags,
		double window_left,
		double window_right,
		int (*fill_event)(
			double trigger_time,
			const std::multimap<double, MapEvent> &match_map,
			FundamentalEvent &fundamental_event,
			std::vector<MatchTriggerStatistics> &statistics
		)
	);

protected:
	// run number
	unsigned int run_;
	// detector name
	std::string name_;

private:
	/// @brief read tirgger from root file
	/// @param[out] trigger_times list of trigger time
	/// @param[in] tag trigger tag to choose trigger file
	/// @returns 0 if success, -1 otherwise
	///
	int ReadTriggerTimes(
		std::vector<double> &trigger_times,
		const std::string &tag = ""
	);


	/// @brief match the map event of detector with trigger
	/// @tparam MapEvent type of detector map event
	/// @param[in] tag tag of trigger to choose file
	/// @param[in] look_window width of look window in nanoseconds
	/// @param[in] hist_look_window histogram to fill time difference
	/// @param[in] window_left left border of matching window in nanoseconds
	/// @param[in] window_right right border of matching window in nanoseconds
	/// @param[out] trigger_times filled trigger times
	/// @param[out] match_map filled match map
	/// @returns entries of detector map tree if success, -1 otherwise
	///
	template<typename MapEvent>
	long long FillMatchMap(
		const std::string &tag,
		double look_window,
		TH1F *hist_look_window,
		double window_left,
		double window_right,
		std::vector<double> &trigger_times,
		std::multimap<double, MapEvent> &match_map
	);
};


template<typename MapEvent>
long long Detector::FillMatchMap(
	const std::string &tag,
	double look_window,
	TH1F *hist_look_window,
	double window_left,
	double window_right,
	std::vector<double> &trigger_times,
	std::multimap<double, MapEvent> &match_map
) {
	// initialize
	trigger_times.clear();
	match_map.clear();

	if (ReadTriggerTimes(trigger_times, tag)) {
		return -1;
	}

	// setup input file and tree
	// name of input file
	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-map-%04u.root",
		kGenerateDataPath, kMappingDir, name_.c_str(), run_
	);
	// pointer to input file
	TFile *ipf = new TFile(input_file_name, "read");
	// pointer to input tree
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// input event
	MapEvent map_event;
	// setup input branches
	map_event.SetupInput(ipt);

	// show begin
	printf("reading %s events   0%%", name_.c_str());
	fflush(stdout);
	// total entries form input tree
	long long entries = ipt->GetEntries();
	// 1/100 of entries for showing process
	long long entry100 = entries / 100 + 1;
	// loop to match trigger and fill events into map
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get entry
		ipt->GetEntry(entry);

		// time that match this event, 0 for not found
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
			// add to histogram for later checks
			hist_look_window->Fill(*iter - map_event.time);
			// jump if out of search window
			if (*iter <= map_event.time + window_left) continue;
			if (*iter >= map_event.time + window_right) continue;

			if (
				fabs(*iter-map_event.time) < fabs(match_time-map_event.time)
			) {
				// found the closer reference trigger time
				match_time = *iter;
			}
		}

		if (match_time != 0.0) {
			// fill to map if found match trigger
			match_map.insert(std::make_pair(match_time, map_event));
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// close input file
	ipf->Close();

	return entries;
}


template<typename MapEvent, typename FundamentalEvent>
int Detector::MatchTrigger(
	const std::string &trigger_tag,
	double window_left,
	double window_right,
	void (*fill_event)(
		double trigger_time,
		const std::multimap<double, MapEvent> &match_map,
		FundamentalEvent &fundamental_event,
		MatchTriggerStatistics &statistics
	)
) {
	const double look_window = 10'000;

	// trigger times
	std::vector<double> trigger_times;
	// histogram of look window
	TH1F *hist_look_window =
		new TH1F("ht", "look window", 1000, -look_window, look_window);
	// detector's match events in map
	std::multimap<double, MapEvent> match_map;
	// detector entries
	long long detector_entries;
	// get detector entries, trigger times and match map
	detector_entries = FillMatchMap<MapEvent>(
		trigger_tag,
		look_window,
		hist_look_window,
		window_left,
		window_right,
		trigger_times,
		match_map
	);
	// return error
	if (detector_entries < 0) return -1;

	// setup output file and tree
	// name of output file
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-fundamental-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		name_.c_str(),
		trigger_tag.empty() ? "" : (trigger_tag + "-").c_str(),
		run_
	);
	// pointer to output file
	TFile *opf = new TFile(output_file_name, "recreate");
	// pointer to output tree
	TTree *opt = new TTree("tree", "tree of fundamental event match trigger");
	// ouput event
	FundamentalEvent fundamental_event;
	// setup output branches
	fundamental_event.SetupOutput(opt);


	// for statistics
	MatchTriggerStatistics statistics(
		run_,
		name_,
		trigger_tag,
		"",
		trigger_times.size(),
		detector_entries
	);

	// show begin
	printf("Writing %s fundamental events   0%%", name_.c_str());
	fflush(stdout);
	// total entries of triggers
	long long entries = trigger_times.size();
	// 1/100 of total triggers
	long long entry100 = entries / 100 + 1;
	// loop for matching trigger time
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// check and fill event
		fill_event(
			trigger_times[entry],
			match_map,
			fundamental_event,
			statistics
		);

		opt->Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// wirte the trees
	opt->Write();
	// write histogram
	hist_look_window->Write();
	// close output file
	opf->Close();

	// save statistics
	statistics.Write();
	// show statistics
	statistics.Print();

	return 0;
}


template<typename FundamentalEvent>
int Detector::VmeMatchTrigger(const std::string &trigger_tag) {
	// read merged vt trigger
	// merged vt fundamental file name
	TString vt_file_name;
	vt_file_name.Form(
		"%s%svt-fundamental-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		trigger_tag.empty() ? "" : (trigger_tag+"-").c_str(),
		run_
	);
	// vt file
	TFile *vt_file = new TFile(vt_file_name, "read");
	// vt tree
	TTree *vt_tree = (TTree*)vt_file->Get("tree");
	if (!vt_tree) {
		std::cerr << "Error: Get tree from "
			<< vt_file_name << " failed.\n";
		return -1;
	}
	// input trigger event
	TriggerEvent vt_event;
	// setup vt tree branches
	vt_event.SetupInput(vt_tree);

	std::vector<long long> vt_times;

	// total number of VME trigger events
	long long entries = vt_tree->GetEntries();
	// 1/100 of entries
	// long long entry100 = entries / 100 + 1;
	// show start
	printf("Reading VME trigger   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		// if (entry % entry100 == 0) {
		// 	printf("\b\b\b\b%3lld%%", entry / entry100);
		// 	fflush(stdout);
		// }
		// get entry
		vt_tree->GetEntry(entry);
		// add valid timestamp or invalid -1 timestamp based on time
		if (vt_event.time > -9e4) vt_times.push_back(vt_event.timestamp);
		else vt_times.push_back(-1);
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// close vt input file
	vt_file->Close();


	// setup detector input and output
	// detector input file name
	TString detector_input_file_name;
	detector_input_file_name.Form(
		"%s%s%s-map-%04u.root",
		kGenerateDataPath, kMappingDir, name_.c_str(), run_
	);
	// detector input file
	TFile *ipf = new TFile(detector_input_file_name, "read");
	// detector input tree
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< detector_input_file_name << " failed.\n";
		return -1;
	}
	// detector output file name
	TString detector_output_file_name;
	detector_output_file_name.Form(
		"%s%s%s-fundamental-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		name_.c_str(),
		trigger_tag.empty() ? "" : (trigger_tag+"-").c_str(),
		run_
	);
	// detector output file
	TFile *opf = new TFile(detector_output_file_name, "recreate");
	// detector output tree
	TTree *opt = new TTree("tree", "fundamental tree");
	// input timestamp
	long long timestamp;
	// input and output event
	FundamentalEvent fundamental_event;
	// setup input and output branches
	ipt->SetBranchAddress("timestamp", &timestamp);
	fundamental_event.SetupInput(ipt);
	fundamental_event.SetupOutput(opt);

	// for statistics
	MatchTriggerStatistics statistics(
		run_,
		name_,
		trigger_tag,
		"",
		vt_times.size(),
		ipt->GetEntries()
	);

	// detecotr map event entries
	long long detector_entries = ipt->GetEntries();
	// detecotr current entry
	long long detector_entry = 0;
	// show start
	printf("Writing %s fundamnetal events   0%%", name_.c_str());
	fflush(stdout);
	for (size_t entry = 0; entry < vt_times.size(); ++entry) {
		// show process
		// if (entry % entry100 == 0) {
		// 	printf("\b\b\b\b%3lld%%", entry / entry100);
		// 	fflush(stdout);
		// }
		
		if (detector_entry < detector_entries) {
			// read detector event
			ipt->GetEntry(detector_entry);
			while (
				timestamp < vt_times[entry]
				&& detector_entry+1 < detector_entries
			) {
				++detector_entry;
				ipt->GetEntry(detector_entry);
			}

			if (timestamp == vt_times[entry]) {
				// record this event and move to next event
				++detector_entry;
				++statistics.used_events;
				++statistics.match_events;
			} else {
				// nullify this event
				fundamental_event.Nullify();
			}
		} else {
			fundamental_event.Nullify();
		}
		// fill event
		opt->Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// write tree
	opt->Write();
	// close files
	opf->Close();
	ipf->Close();

	// save statistics
	statistics.Write();
	// show statistics
	statistics.Print();

	return 0;
}


template<typename MapEvent, typename FundamentalEvent>
int Detector::ExtractTrigger(
	const std::string &trigger_tag,
	const std::vector<std::string> &extract_tags,
	double window_left,
	double window_right,
	int (*fill_event)(
		double trigger_time,
		const std::multimap<double, MapEvent> &match_map,
		FundamentalEvent &fundamental_event,
		std::vector<MatchTriggerStatistics> &statistics
	)
) {
	const double look_window = 10'000;

	// trigger times
	std::vector<double> trigger_times;
	// histogram for look window
	TH1F *hist_look_window =
		new TH1F("ht", "look window", 1000, -look_window, look_window);
	// detector's match events in map
	std::multimap<double, MapEvent> match_map;
	// detector entries
	long long detector_entries;
	// get detector entries, trigger times and match map
	detector_entries = FillMatchMap<MapEvent>(
		trigger_tag,
		look_window,
		hist_look_window,
		window_left,
		window_right,
		trigger_times,
		match_map
	);
	// return error
	if (detector_entries < 0) return -1;

	// output files and trees
	// output extracted trigger files
	std::vector<TFile*> out_trigger_files;
	// output detector fundamental event files
	std::vector<TFile*> out_detector_files;
	// output extracted trigger trees
	std::vector<TTree*> out_trigger_trees;
	// output detector fundamental event trees
	std::vector<TTree*> out_detector_trees;

	// output data
	// output trigger event
	double xia_trigger_time;
	// output detector event
	FundamentalEvent fundamental_event;


	for (const std::string &tag : extract_tags) {
		// setup output extracted trigger files and trees
		// output extracted trigger file name
		TString out_trigger_file_name;
		out_trigger_file_name.Form(
			"%s%sxt-map-%s%s-%04u.root",
			kGenerateDataPath,
			kMappingDir,
			trigger_tag.empty() ? "" : (trigger_tag + "-").c_str(),
			tag.c_str(),
			run_
		);
		out_trigger_files.push_back(
			new TFile(out_trigger_file_name, "recreate")
		);
		// output extracted trigger tree
		TTree* out_trigger_tree =
			new TTree("tree", "tree of extracted trigger");
		out_trigger_trees.push_back(out_trigger_tree);
		// setup output trigger branch
		out_trigger_tree->Branch("time", &xia_trigger_time, "t/D");

		// setup output detector fundamental file and tree
		// output detector file name
		TString out_detector_file_name;
		out_detector_file_name.Form(
			"%s%s%s-fundamental-%s%s-%04u.root",
			kGenerateDataPath,
			kFundamentalDir,
			name_.c_str(),
			trigger_tag.empty() ? "" : (trigger_tag + "-").c_str(),
			tag.c_str(),
			run_
		);
		// output fundamental file
		out_detector_files.push_back(
			new TFile(out_detector_file_name, "recreate")
		);
		// output detector tree
		TTree *out_detector_tree =
			new TTree("tree", "tree of fundamental event match trigger");
		out_detector_trees.push_back(out_detector_tree);
		// setup output detector branch
		fundamental_event.SetupOutput(out_detector_tree);
	}

	// statistics
	std::vector<MatchTriggerStatistics> statistics;
	for (const std::string &tag : extract_tags) {
		statistics.emplace_back(
			run_,
			name_,
			trigger_tag,
			tag,
			trigger_times.size(),
			detector_entries
		);
	}

	// show begin
	printf("Writing %s fundamental events   0%%", name_.c_str());
	fflush(stdout);
	// total number of triggers
	long long entries = trigger_times.size();
	// 1/100 of total entries
	long long entry100 = entries / 100 + 1;
	// loop for matching trigger time
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// check and fill event
		int fill_index = fill_event(
			trigger_times[entry],
			match_map,
			fundamental_event,
			statistics
		);

		xia_trigger_time = trigger_times[entry];
		if (fill_index >= 0) {
			out_trigger_trees[fill_index]->Fill();
			out_detector_trees[fill_index]->Fill();
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	for (size_t i = 0; i < out_trigger_files.size(); ++i) {
		// write trigger tree and close file
		out_trigger_files[i]->cd();
		out_trigger_trees[i]->Write();
		out_trigger_files[i]->Close();
		// wirte detector tree, look window and close file
		out_detector_files[i]->cd();
		out_detector_trees[i]->Write();
		hist_look_window->Write();
		out_detector_files[i]->Close();
	}

	for (MatchTriggerStatistics &s : statistics) {
		// save statistics
		s.Write();
		// show statistics
		s.Print();
	}

	return 0;
}


// //-----------------------------------------------------------------------------
// //									merge
// //-----------------------------------------------------------------------------
// public:
// 	/// @brief merge events in adjacent strip
// 	///
// 	/// @param[in] energy_cut front back energy cut
// 	/// @returns 0 if success, -1 otherwise
// 	///
// 	virtual int Merge(double energy_cut = 0.02);

// 	// correlated and merged event
// 	struct MergedEvent
// 	{
// 		unsigned short hit;
// 		long long timestamp;
// 		double front_strip[4];
// 		double back_strip[4];
// 		double time[4];
// 		double energy[4];
// 	};

// protected:

// 	/// @brief setup merged output file and tree
// 	///
// 	/// @param[out] output_file pointer to pointer to the output file
// 	/// @param[out] output_tree pointer to pointer to the output tree
// 	///
// 	virtual void SetupMergedOutput(TFile **output_file, TTree **output_tree);


// 	/// @brief setup unmerged output file
// 	///
// 	/// @param[out] output_file pointer to pointer to unmerged file
// 	/// @param[out] output_tree pointer to pointer to unmerged tree
// 	///
// 	virtual void SetupUnmergedOutput(
// 		TFile **output_file,
// 		TTree **output_tree
// 	);


// 	// merged event
// 	MergedEvent merged_;


// //-----------------------------------------------------------------------------
// //									telescope
// //-----------------------------------------------------------------------------
// public:

// 	/// @brief read merged events
// 	///
// 	/// @returns 0 if success, -1 otherwise
// 	///
// 	virtual int ReadMergedEvents();


// 	/// @brief get total merged event number
// 	///
// 	/// @returns merged event number
// 	///
// 	virtual inline long long GetEvents()
// 	{
// 		return merged_events_;
// 	}


// 	/// @brief get pointer to merged event by index
// 	///
// 	/// @param[in] index index of event to get
// 	/// @returns pointer to merged-event, or nullptr if not found
// 	///
// 	virtual inline MergedEvent* GetEvent(long long index) {
// 		if (index >= merged_events_ || index < 0)
// 		{
// 			return nullptr;
// 		}
// 		input_merged_tree_->GetEntry(index);
// 		return &merged_;
// 	}


// 	/// @brief finish reading merged-events
// 	virtual inline void FinishMergedEvents()
// 	{
// 		input_merged_file_->Close();
// 	}

// private:
// 	// pointer to input merged file
// 	TFile *input_merged_file_;
// 	// pointer to input merged tree
// 	TTree *input_merged_tree_;
// 	// total number of merged-events
// 	long long merged_events_;



}

#endif			// __DETECTOR_H__