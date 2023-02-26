#ifndef __DETECTOR_H__
#define __DETECTOR_H__

#include <vector>
#include <iostream>

#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

#include "include/defs.h"

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
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @returns 0 if success, -1 otherwise
	///
	virtual int MatchTrigger(double window_left, double window_right) = 0;


	/// @brief template of match trigger
	/// @tparam MapEvent type of map event
	/// @tparam FundamentalEvent type of fundamental event
	/// @tparam Statistics type of statistics
	/// @param[in] window_left left edge of match window
	/// @param[in] window_right right edge of match window
	/// @param[in] fill_event function to convert map event
	///		to fundamental event
	/// @returns 0 if success, -1 otherwise
	///
	template<
		typename MapEvent,
		typename FundamentalEvent,
		typename Statistics
	>
	int MatchTrigger(
		double window_left,
		double window_right,
		void (*fill_event)(
			double trigger_time,
			const std::multimap<double, MapEvent> &match_map,
			FundamentalEvent &fundamental_event,
			Statistics &statistics
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
	/// @returns 0 if success, -1 otherwise
	///
	int ReadTriggerTimes(std::vector<double> &trigger_times);
};



template<
	typename MapEvent,
	typename FundamentalEvent,
	typename Statistics
>
int Detector::MatchTrigger(
	double window_left,
	double window_right,
	void (*fill_event)(
		double trigger_time,
		const std::multimap<double, MapEvent> &match_map,
		FundamentalEvent &fundamental_event,
		Statistics &statistics
	)
) {
	const double look_window = 10'000;

	std::vector<double> trigger_times;
	if (ReadTriggerTimes(trigger_times)) {
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

	// setup output file and tree
	// name of output file
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-fundamental-%04u.root",
		kGenerateDataPath, kFundamentalDir, name_.c_str(), run_
	);
	// pointer to output file
	TFile *opf = new TFile(output_file_name, "recreate");
	// histogram of look window
	TH1F *hist_look_window =
		new TH1F("ht", "look window", 1000, -look_window, look_window);
	// pointer to output tree
	TTree *opt = new TTree("tree", "tree of fundamental event match trigger");
	// ouput event
	FundamentalEvent fundamental_event;
	// setup output branches
	fundamental_event.SetupOutput(opt);

	// detector's match events in map
	std::multimap<double, MapEvent> match_map;

	// show begin
	printf("reading detector events   0%%");
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
	// write histogram
	hist_look_window->Write();
	// close input file
	ipf->Close();

	// for statistics
	Statistics statistics(trigger_times.size());

	// show begin
	printf("Writing fundamental events   0%%");
	fflush(stdout);
	entries = trigger_times.size();
	entry100 = entries / 100 + 1;
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
	// close output file
	opf->Close();

	// show statistics
	std::cout << statistics << "\n";

	return 0;
}

// //-----------------------------------------------------------------------------
// //										geometry
// //-----------------------------------------------------------------------------
// protected:
// 	/// @brief get front strip strip
// 	///
// 	/// @returns front strip number
// 	///
// 	virtual inline size_t FrontStrip()
// 	{
// 		return 0;
// 	}


// 	/// @brief get back strip number
// 	///
// 	/// @returns back strip number
// 	///
// 	virtual inline size_t BackStrip()
// 	{
// 		return 0;
// 	}



// //-----------------------------------------------------------------------------
// //						single side and correlated events
// //-----------------------------------------------------------------------------
// public:

// 	/// @brief correlate events in detector
// 	///
// 	///	@param[in] single_side_window single side time window
// 	/// @param[in] double_sides_window double side time window
// 	/// @returns 0 if success, -1 otherwise
// 	///
// 	virtual int Correlate
// 	(
// 		unsigned int single_side_window,
// 		unsigned int double_sides_window
// 	);

// protected:
// 	// correlated event
// 	struct CorrelatedEvent
// 	{
// 		unsigned short index;
// 		unsigned short front_hit;
// 		unsigned short back_hit;
// 		long long timestamp;
// 		unsigned short front_strip[8];
// 		unsigned short back_strip[8];
// 		double front_time[8];
// 		double back_time[8];
// 		double front_energy[8];
// 		double back_energy[8];
// 	};


// 	/// @brief setup residual single side output file and tree
// 	///
// 	/// @param[out] output_file pointer to pointer to output file
// 	/// @param[out] output_tree pointer to pointer to output tree
// 	///
// 	virtual void SetupSingleSideOutput
// 	(
// 		TFile **output_file,
// 		TTree **output_tree
// 	);


// 	/// @brief setup residual correlated output file and tree
// 	///
// 	/// @param[out] output_file pointer to pointer to correlated file
// 	/// @param[out] output_tree pointer to pointer to correalted tree
// 	///
// 	virtual void SetupCorrelatedOutput(
// 		TFile **output_file,
// 		TTree **output_tree
// 	);


// 	/// @brief abstrct function setup output tree of correlated events
// 	///
// 	/// @param[in] tree output correlated tree to setup
// 	///
// 	virtual void SetupOutputCorrelatedTree(TTree *tree) = 0;


// 	/// @brief read correlated events from root file
// 	///
// 	/// @returns 0 for success, -1 otherwise
// 	///
// 	virtual void SetupInputCorrelatedTree(TTree *tree) = 0;


// 	// residual correlated event
// 	CorrelatedEvent correlation_;

// private:
// 	// single event before correlation
// 	struct Event
// 	{
// 		bool used;
// 		unsigned short index;
// 		unsigned short side;
// 		unsigned short strip;
// 		double time;
// 		double energy;
// 	};


// 	// single side event after single-side-correlation before double-sides-correlation
// 	struct SingleSideEvent
// 	{
// 		bool correlated;
// 		unsigned short index;
// 		unsigned short side;
// 		unsigned short hit;
// 		long long timestamp;
// 		unsigned short strip[8];
// 		double time[8];
// 		double energy[8];
// 	};


// 	// information from read events
// 	struct SingleSideInfo
// 	{
// 		TFile *single_side_file;
// 		TTree *single_side_tree;
// 		std::multimap<long long, SingleSideEvent> events;
// 	};


// 	/// @brief read detector events from mapped file and store in map
// 	///
// 	/// @param[in] single_side_window single side time window
// 	/// @param[out] info output single event information
// 	/// @returns 0 for success, -1 otherwise
// 	///
// 	int ReadEvents(unsigned int single_side_window, SingleSideInfo &info);


// 	/// @brief store residual single-side events
// 	///
// 	/// @param[in] info information about single side events
// 	///
// 	void StoreResidualSingleSideEvents(const SingleSideInfo &info);


// 	/// @brief sort single-side-events by strip
// 	///
// 	/// @param[in] event pointer to single-side-event
// 	///
// 	void BubbleSort(SingleSideEvent *event);


// 	// single side event
// 	SingleSideEvent single_side_;


// //-----------------------------------------------------------------------------
// //									normalize
// //-----------------------------------------------------------------------------
// public:
// 	/// @brief calculate normalize parameters
// 	///
// 	/// @param[in] length number of run to chain
// 	/// @param[in] ref_front reference front strip
// 	/// @param[in] ref_back reference back strip
// 	/// @param[in] iteration run in iteration mode
// 	/// @returns 0 for succes, -1 otherwise
// 	///
// 	virtual int Normalize
// 	(
// 		int length,
// 		unsigned short ref_front,
// 		unsigned short ref_back,
// 		bool iteration
// 	);

// 	/// @brief build normalize result
// 	///
// 	virtual int NormalCorrelate();

// protected:

// 	/// @brief calculate normalize parameters of first side
// 	///
// 	/// @param[in] chain input tree
// 	/// @param[in] ref_front reference front strip
// 	/// @param[in] ref_back reference back strip
// 	/// @param[in] g_front_back_energy list of graph of front-back energy
// 	/// @param[in] g_back_front_energy list of graph of back-front energy
// 	/// @param[in] iteration iteration mode
// 	///
// 	virtual void NormalizeFirstSide
// 	(
// 		TTree *chain,
// 		unsigned short ref_front,
// 		unsigned short ref_back,
// 		TGraph **g_front_back_energy,
// 		TGraph **g_back_front_energy,
// 		bool iteration
// 	);


// 	/// @brief calculate normalize parameters of first side
// 	///
// 	/// @param[in] chain input tree
// 	/// @param[in] ref_front reference front strip
// 	/// @param[in] ref_back reference back strip
// 	/// @param[in] g_front_back_energy list of graph of front-back energy
// 	/// @param[in] g_back_front_energy list of graph of back-front energy
// 	/// @param[in] iteration iteration mode
// 	///
// 	virtual void NormalizeSecondSide
// 	(
// 		TTree *chain,
// 		unsigned short ref_front,
// 		unsigned short ref_back,
// 		TGraph **g_front_back_energy,
// 		TGraph **g_back_front_energy,
// 		bool iteration
// 	);


// 	/// @brief check whether the energy of correlated event can be use in
// 	/// 	normalization for back strip
// 	///
// 	/// @param[in] correlation corrlated event
// 	/// @param[in] iteration iteration mode
// 	/// @returns true if pass, or false otherwise
// 	///
// 	virtual bool NormalizeFrontEnergyCheck
// 	(
// 		const CorrelatedEvent &correlation,
// 		bool iteration
// 	) = 0;


// 	/// @brief check whether the energy of correlated event can be use in
// 	///		normalization for front strip
// 	///
// 	/// @param[in] correlation correlated event
// 	/// @param[in] iteration iteration mode
// 	/// @returns ture if pass, or false otherwise
// 	///
// 	virtual bool NormalizeBackEnergyCheck
// 	(
// 		const CorrelatedEvent &correlation,
// 		bool iteration
// 	) = 0;


// 	/// @brief calculate normalized energy
// 	///
// 	/// @param[in] side 0 for front, 1 for back
// 	/// @param[in] strip strip number
// 	/// @param[in] energy energy before normalization
// 	/// @returns energy after normalization
// 	///
// 	inline double NormalizeEnergy
// 	(
// 		size_t side,
// 		unsigned short strip,
// 		double energy
// 	)
// 	{
// 		if (side == 0)
// 		{
// 			return normalize_front_param_[strip][0]
// 				+ normalize_front_param_[strip][1] * energy;
// 		}
// 		else
// 		{
// 			return normalize_back_param_[strip][0]
// 				+ normalize_back_param_[strip][1] * energy;
// 		}
// 	}

// 	// normalize back parameters
// 	double normalize_back_param_[64][2];
// 	// normalize front parameters
// 	double normalize_front_param_[64][2];

// private:
// 	/// @brief read normalize parameters from file
// 	///
// 	/// @returns 0 for success, -1 otherwise
// 	///
// 	int ReadNormalizeParameters();


// 	/// @brief write normalize parameters to file
// 	///
// 	/// @returns 0 for success, -1 otherwise
// 	///
// 	int WriteNormalizeParameters();

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