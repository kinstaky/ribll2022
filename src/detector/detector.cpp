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
	const std::string &name
)
: run_(run)
, name_(name)
{

	// for (size_t i = 0; i < 64; ++i)
	// {
	// 	normalize_front_param_[i][0] = 0.0;
	// 	normalize_front_param_[i][1] = 1.0;
	// 	normalize_back_param_[i][0] = 0.0;
	// 	normalize_back_param_[i][1] = 1.0;
	// }
}

int Detector::ReadTriggerTimes(std::vector<double> &trigger_times) {
	// clear data
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

Detector::MatchTriggerStatistics::MatchTriggerStatistics(
	long long reference,
	long long total
)
: reference_events(reference)
, total_events(total)
, match_events(0)
, used_events(0)
, oversize_events(0) {
}


std::ostream& operator<<(
	std::ostream &os,
	const Detector::MatchTriggerStatistics &st
) {
	return os << "trigger match rate "
		<< st.match_events << " / " << st.reference_events << "  "
		<< double(st.match_events) / double(st.reference_events) << "\n"
		<< "mapped event used rate "
		<< st.used_events << " / " << st.total_events << "  "
		<< double(st.used_events) / double(st.total_events) << "\n"
		<< "oversize rate "
		<< st.oversize_events << " / " << st.reference_events << "  "
		<< double(st.oversize_events) / double(st.reference_events);
}


// //-----------------------------------------------------------------------------
// //						single side and correlated events
// //-----------------------------------------------------------------------------

// void Detector::SetupSingleSideOutput(TFile **output_file, TTree **output_tree)
// {
// 	// single side event file name
// 	TString single_side_file_name;
// 	single_side_file_name.Form(
// 		"%s%s%s-single-side-%04d.root",
// 		kGenerateDataPath, kSingleSideDir, name_.c_str(), run_
// 	);
// 	// setup single side file
// 	TFile *single_side_file = new TFile(single_side_file_name, "recreate");
// 	// setup single side tree
// 	TTree *single_side_tree = new TTree(
// 		"tree",
// 		TString::Format("residual single side events of %s", name_.c_str())
// 	);
// 	// setup single side tree branches
// 	single_side_tree->Branch("timestamp", &single_side_.timestamp, "ts/L");
// 	single_side_tree->Branch("index", &single_side_.index, "index/s");
// 	single_side_tree->Branch("side", &single_side_.side, "side/s");
// 	single_side_tree->Branch("hit", &single_side_.hit, "hit/s");
// 	single_side_tree->Branch("strip", single_side_.strip, "s[hit]/s");
// 	single_side_tree->Branch("time", single_side_.time, "t[hit]/D");
// 	single_side_tree->Branch("energy", single_side_.energy, "e[hit]/D");

// 	// output pointer to single side file and single side tree
// 	*output_file = single_side_file;
// 	*output_tree = single_side_tree;

// 	return;
// }


// int Detector::ReadEvents(unsigned int single_side_window, SingleSideInfo &info)
// {
// 	// setup input file
// 	TString input_file_name;
// 	input_file_name.Form(
// 		"%s%s%s-map-%04d.root",
// 		kGenerateDataPath,
// 		kMappingDir,
// 		name_.c_str(),
// 		run_
// 	);
// 	TFile *input_file = new TFile(input_file_name, "read");
// 	TTree *ipt = (TTree*)input_file->Get("tree");
// 	if (!ipt)
// 	{
// 		std::cerr << "Error: get tree from " << input_file_name
// 			<< " failed.\n";
// 		return -1;
// 	}
// 	// input data and branches
// 	long long timestamp;
// 	unsigned short index, side, strip;
// 	double time, energy;
// 	// branches
// 	ipt->SetBranchAddress("timestamp", &timestamp);
// 	ipt->SetBranchAddress("index", &index);
// 	ipt->SetBranchAddress("side", &side);
// 	ipt->SetBranchAddress("strip", &strip);
// 	ipt->SetBranchAddress("time", &time);
// 	ipt->SetBranchAddress("energy", &energy);

// 	SetupSingleSideOutput(&info.single_side_file, &info.single_side_tree);
// 	TH1F *single_side_look_window = new TH1F(
// 		"ht", "look window", 1000, -look_window, look_window
// 	);
// 	TH1F *front_side_look_window = new TH1F(
// 		"hft", "look window", 1000, -look_window, look_window
// 	);
// 	TH1F *back_side_look_window = new TH1F(
// 		"hbt", "look window", 1000, -look_window, look_window
// 	);

// 	// initialize
// 	info.events.clear();

// 	// show process
// 	printf("Reading events   0%%");
// 	fflush(stdout);
// 	// entries
// 	long long entries = ipt->GetEntries();
// 	// 1/100 of entry
// 	long long entry100 = entries / 100;

// 	// to show hit rate at the end of this function
// 	// total event number of front side
// 	long long front_total_event = 0;
// 	// total event number of back side
// 	long long back_total_event = 0;
// 	// front side single event counts
// 	long long front_single_event = 0;
// 	// back side single event counts
// 	long long back_single_event = 0;
// 	// front side double event counts
// 	long long front_double_event = 0;
// 	// back side double event counts
// 	long long back_double_event = 0;
// 	// front side triple event counts
// 	long long front_triple_event = 0;
// 	// back side triple event counts
// 	long long back_triple_event = 0;

// 	for (long long entry = 0; entry < entries; ++entry)
// 	{
// 		if (entry % entry100 == 0)
// 		{
// 			printf("\b\b\b\b%3lld%%", entry / entry100);
// 			fflush(stdout);
// 		}

// 		ipt->GetEntry(entry);
// 		if (energy < 100) continue;

// 		// iterator points to the element to erase, to make the key be the
// 		// timestmap of the event with largest energy
// 		auto erase_iter = info.events.end();
// 		// search if in the time-window of exist event
// 		bool found = false;
// 		for
// 		(
// 			auto iter = info.events.lower_bound(timestamp - look_window);
// 			iter != info.events.lower_bound(timestamp + look_window);
// 			++iter
// 		)
// 		{
// 			// ignore other detector event (different index)
// 			if (iter->second.index != index) continue;
// 			// ignore other side event
// 			if (iter->second.side != side) continue;
// 			// fill to look window
// 			single_side_look_window->Fill(iter->first - timestamp);
// 			if (!side)
// 			{
// 				front_side_look_window->Fill(iter->first - timestamp);
// 			}
// 			else
// 			{
// 				back_side_look_window->Fill(iter->first - timestamp);
// 			}
// 			// continue if out of single side window
// 			if (iter->first <= timestamp - single_side_window) continue;
// 			if (iter->first >= timestamp + single_side_window) continue;
// 			// jump if has found correlated event
// 			if (found) continue;
// 			// found correlated event
// 			found = true;
// 			if (iter->second.hit < 8)
// 			{
// 				iter->second.strip[iter->second.hit] = strip;
// 				iter->second.time[iter->second.hit] = time;
// 				iter->second.energy[iter->second.hit] = energy;
// 				// use the largest energy event's timestamp
// 				bool energy_max = true;
// 				for (size_t i = 0; i < iter->second.hit; ++i)
// 				{
// 					if (energy < iter->second.energy[i])
// 					{
// 						energy_max = false;
// 						break;
// 					}
// 				}
// 				// new event has the largest erergy
// 				// use its timestamp to represent the whole event
// 				if (energy_max)
// 				{
// 					erase_iter = iter;
// 				}
// 			}
// 			++(iter->second.hit);
// 			// record single, double and triple event counts
// 			switch (iter->second.hit)
// 			{
// 				case 2:
// 					if (iter->second.side == 0)
// 					{
// 						// front side
// 						--front_single_event;
// 						++front_double_event;
// 					}
// 					else
// 					{
// 						// back side
// 						--back_single_event;
// 						++back_double_event;
// 					}
// 					break;
// 				case 3:
// 					if (iter->second.side == 0)
// 					{
// 						// front side
// 						--front_double_event;
// 						++front_triple_event;
// 					}
// 					else
// 					{
// 						// back side
// 						--back_double_event;
// 						++back_triple_event;
// 					}
// 					break;
// 				default:
// 					// over 3
// 					if (iter->second.side == 0)
// 					{
// 						// front side
// 						--front_triple_event;
// 					}
// 					else
// 					{
// 						// back side
// 						--back_triple_event;
// 					}
// 			}
// 		}
// 		if (found && erase_iter != info.events.end())
// 		{
// 			SingleSideEvent event = erase_iter->second;
// 			event.timestamp = timestamp;
// 			info.events.erase(erase_iter);
// 			info.events.insert(std::make_pair(timestamp, event));
// 		}
// 		if (!found)
// 		{
// 			// not found, insert new element to map
// 			SingleSideEvent event{
// 				false,
// 				index,
// 				side,
// 				1,
// 				timestamp,
// 				{strip, 0, 0, 0, 0, 0, 0, 0},
// 				{time, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
// 				{energy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
// 			};
// 			info.events.insert(std::make_pair(timestamp, event));
// 			// increase single event
// 			if (side == 0)
// 			{
// 				// front side
// 				++front_single_event;
// 				++front_total_event;
// 			}
// 			else
// 			{
// 				// back side
// 				++back_single_event;
// 				++back_total_event;
// 			}
// 		}
// 	}
// 	printf("\b\b\b\b100%%\n");

// 	// write look window
// 	info.single_side_file->cd();
// 	single_side_look_window->Write();
// 	front_side_look_window->Write();
// 	back_side_look_window->Write();
// 	// close input file
// 	input_file->Close();

// 	// show rate
// 	std::cout
// 		<< "front side total " << front_total_event << "\n"
// 		<< "  single rate " << front_single_event << " "
// 		<< double(front_single_event) / front_total_event << "\n"
// 		<< "  double rate " << front_double_event << " "
// 		<< double(front_double_event) / front_total_event << "\n"
// 		<< "  triple rate " << front_triple_event << " "
// 		<< double(front_triple_event) / front_total_event << "\n"
// 		<< "back side total " << back_total_event << "\n"
// 		<< "  single rate " << back_single_event << " "
// 		<< double(back_single_event) / back_total_event << "\n"
// 		<< "  double rate " << back_double_event << " "
// 		<< double(back_double_event) / back_total_event << "\n"
// 		<< "  triple rate " << back_triple_event << " "
// 		<< double(back_triple_event) / back_total_event << "\n";

// 	return 0;
// }


// void Detector::BubbleSort(SingleSideEvent *event)
// {
// 	bool exchange = true;
// 	while (exchange)
// 	{
// 		exchange = false;
// 		for (unsigned short i = 0; i < event->hit-1; ++i)
// 		{
// 			if (event->strip[i] > event->strip[i+1])
// 			{
// 				// exchange
// 				exchange = true;
// 				// exchange strip
// 				unsigned short tmp_strip = event->strip[i];
// 				event->strip[i] = event->strip[i+1];
// 				event->strip[i+1] = tmp_strip;
// 				// excahange time
// 				double tmp = event->time[i];
// 				event->time[i] = event->time[i+1];
// 				event->time[i+1] = tmp;
// 				// exchange energy
// 				tmp = event->energy[i];
// 				event->energy[i] = event->energy[i+1];
// 				event->energy[i+1] = tmp;
// 			}
// 		}
// 	}
// 	return;
// }


// void Detector::StoreResidualSingleSideEvents(const SingleSideInfo &info) {
// 	// fill residual single side events
// 	// show process
// 	printf("\nFilling residual single side events   0%%");
// 	fflush(stdout);
// 	// filtering index
// 	size_t filtering_index = 0;
// 	// 1/100 of events size
// 	size_t filtering_index_100 = info.events.size() / 100;
// 	info.single_side_file->cd();
// 	for (const auto &event : info.events)
// 	{
// 		// show process
// 		if (filtering_index % filtering_index_100 == 0) {
// 			printf("\b\b\b\b%3lu%%", filtering_index / filtering_index_100);
// 			fflush(stdout);
// 		}
// 		++filtering_index;

// 		if (event.second.correlated) continue;

// 		single_side_ = event.second;

// 		if (single_side_.hit > 8) continue;
// 		// bubble sort by strip
// 		BubbleSort(&single_side_);

// 		info.single_side_tree->Fill();
// 	}
// 	// show process
// 	printf("\b\b\b\b100%%\n");
// 	// write tree
// 	info.single_side_tree->Write();
// }


// void Detector::SetupCorrelatedOutput
// (
// 	TFile **output_file,
// 	TTree **output_tree
// )
// {
// 	// setup output correlated file
// 	TString correlated_file_name;
// 	correlated_file_name.Form(
// 		"%s%s%s-correlated-%04d.root",
// 		kGenerateDataPath, kCorrelationDir, name_.c_str(), run_
// 	);
// 	*output_file = new TFile(correlated_file_name, "recreate");
// 	// setup output tree
// 	*output_tree = new TTree(
// 		"tree", TString::Format("correlated tree of %s", name_.c_str())
// 	);
// 	// setup output data
// 	SetupOutputCorrelatedTree(*output_tree);

// 	return;
// }



// int Detector::Correlate
// (
// 	unsigned int single_side_window,
// 	unsigned int double_sides_window
// )
// {
// 	SingleSideInfo info;
// 	info.single_side_file = nullptr;
// 	info.single_side_tree = nullptr;
// 	// read detector events
// 	if (ReadEvents(single_side_window, info))
// 	{
// 		std::cerr << "Error: Read " << name_ << " events failed.\n";
// 		return -1;
// 	}


// 	TFile *correlated_file = nullptr;
// 	TTree *correlated_tree = nullptr;
// 	SetupCorrelatedOutput(&correlated_file, &correlated_tree);
// 	correlated_file->cd();
// 	TH1F *double_side_look_window
// 		= new TH1F("ht", "look window", 1000, 0, look_window);

// 	// show process
// 	printf("\nCorrelating %s   0%%", name_.c_str());
// 	fflush(stdout);
// 	// 1/100 of entry
// 	long long entry100 = info.events.size() / 100;
// 	// correlating entry
// 	long long entry = 0;

// 	// for calculating the correlated rate
// 	// correlated rate of events
// 	long long correlated_events = 0;
// 	// over correlated events, one to serveral condition
// 	long long over_correlated_events = 0;
// 	// used events
// 	long long used_events = 0;

// 	for
// 	(
// 		auto ievent = info.events.begin();
// 		ievent != info.events.end();
// 		++ievent
// 	)
// 	{
// 		// show process
// 		if (entry % entry100 == 0)
// 		{
// 			printf("\b\b\b\b%3lld%%", entry / entry100);
// 			fflush(stdout);
// 		}
// 		++entry;

// 		// jump correlated events
// 		if (ievent->second.correlated) continue;
// 		// jump if hit over 8
// 		if (ievent->second.hit >= 8) continue;

// 		// found correlated event on the other side
// 		bool correlated = false;
// 		for
// 		(
// 			auto jevent = info.events.lower_bound(ievent->first);
// 			jevent != info.events.upper_bound(ievent->first + look_window);
// 			++jevent
// 		)
// 		{
// 			// ignore itself
// 			if (jevent == ievent) continue;
// 			// fill to look window
// 			double_side_look_window->Fill(jevent->first - ievent->first);

// 			// jump if out of search window
// 			if (jevent->first - ievent->first > double_sides_window) continue;

// 			// jump other detector
// 			if (jevent->second.index != ievent->second.index) continue;
// 			// jump same side event
// 			if (jevent->second.side == ievent->second.side) continue;
// 			// jump if hit over 8
// 			if (jevent->second.hit >= 8) continue;

// 			// jump correlated event
// 			if (jevent->second.correlated)
// 			{
// 				++used_events;
// 				continue;
// 			}

// 			if (!correlated)
// 			{
// 				// x side event
// 				SingleSideEvent *front_event = nullptr;
// 				// y side event
// 				SingleSideEvent *back_event = nullptr;
// 				if (ievent->second.side == 0)
// 				{
// 					front_event = &(ievent->second);
// 					back_event = &(jevent->second);
// 					correlation_.timestamp = ievent->first;
// 				}
// 				else
// 				{
// 					front_event = &(jevent->second);
// 					back_event = &(ievent->second);
// 					correlation_.timestamp = jevent->first;
// 				}

// 				BubbleSort(front_event);
// 				BubbleSort(back_event);

// 				correlated = true;
// 				jevent->second.correlated = true;

// 				// fill to correlation_
// 				correlation_.index = front_event->index;
// 				// sort event in single side by energy
// 				correlation_.front_hit = front_event->hit;
// 				for (unsigned short i = 0; i < front_event->hit; ++i)
// 				{
// 					correlation_.front_strip[i] = front_event->strip[i];
// 					correlation_.front_time[i] = front_event->time[i];
// 					correlation_.front_energy[i] = front_event->energy[i];
// 				}
// 				correlation_.back_hit = back_event->hit;
// 				for (unsigned short i = 0; i < back_event->hit; ++i)
// 				{
// 					correlation_.back_strip[i] = back_event->strip[i];
// 					correlation_.back_time[i] = back_event->time[i];
// 					correlation_.back_energy[i] = back_event->energy[i];
// 				}
// 			}
// 		}

// 		if (correlated)
// 		{
// 			ievent->second.correlated = true;
// 			correlated_tree->Fill();
// 			correlated_events += 2;
// 		}

// 	}
// 	// show finished
// 	printf("\b\b\b\b100%%\n");

// 	// correlated events
// 	correlated_file->cd();
// 	double_side_look_window->Write();
// 	correlated_tree->Write();
// 	correlated_file->Close();

// 	// show correlated rate
// 	std::cout
// 		<< "correlated rate " << correlated_events << " / "
// 		<< info.events.size() << " "
// 		<< double(correlated_events) / double(info.events.size()) << "\n"
// 		<< "over correlated " << over_correlated_events << "\n"
// 		<< "used events " << used_events << "\n";


// 	StoreResidualSingleSideEvents(info);

// 	return 0;
// }




// //-----------------------------------------------------------------------------
// //									normalize
// //-----------------------------------------------------------------------------

// int Detector::ReadNormalizeParameters()
// {
// 	// setup file name
// 	TString file_name;
// 	file_name.Form
// 	(
// 		"%s%s%s.txt",
// 		kGenerateDataPath, kNormalizeDir, name_.c_str()
// 	);
// 	// open file
// 	std::ifstream fin(file_name.Data());
// 	if (!fin.good())
// 	{
// 		std::cerr << "Error: open normalize file "
// 			<< file_name << " failed.\n";
// 		return -1;
// 	}

// 	// setup strip number
// 	size_t strip_num;

// 	// read front strip number and normalize parameters
// 	fin >> strip_num;
// 	// read normalized paramters for front strips
// 	for (size_t i = 0; i < strip_num; ++i)
// 	{
// 		fin >> normalize_front_param_[i][0] >> normalize_front_param_[i][1];
// 	}

// 	//  read back strip number
// 	fin >> strip_num;
// 	// read normalized parameters for back strips
// 	for (size_t i = 0; i < strip_num; ++i)
// 	{
// 		fin >> normalize_back_param_[i][0] >> normalize_back_param_[i][1];
// 	}

// 	// close file
// 	fin.close();

// 	return 0;
// }


// int Detector::WriteNormalizeParameters()
// {
// 	// setup file name
// 	TString file_name;
// 	file_name.Form
// 	(
// 		"%s%s%s.txt",
// 		kGenerateDataPath, kNormalizeDir, name_.c_str()
// 	);
// 	// open output file
// 	std::ofstream fout(file_name.Data());
// 	if (!fout.good())
// 	{
// 		std::cerr << "Error: open normalize parameters file "
// 			<< file_name << " failed.\n";
// 		return -1;
// 	}

// 	// write front strip number
// 	fout << FrontStrip() << "\n";
// 	// write normalized paramters for front strips
// 	for (size_t i = 0; i < FrontStrip(); ++i)
// 	{
// 		fout << normalize_front_param_[i][0]<< " "
// 			<< normalize_front_param_[i][1] << "\n";
// 	}
// 	// write back strip number
// 	fout << BackStrip() << "\n";
// 	// write normalized parameters for back strips
// 	for (size_t i = 0; i < BackStrip(); ++i)
// 	{
// 		fout << normalize_back_param_[i][0] << " "
// 			<< normalize_back_param_[i][1] << "\n";
// 	}

// 	// close file
// 	fout.close();

// 	return 0;
// }


// void Detector::NormalizeFirstSide(
// 	TTree *chain,
// 	unsigned short,
// 	unsigned short ref_back,
// 	TGraph **g_front_back_energy,
// 	TGraph **,
// 	bool iteration
// ) {
// 	// show process
// 	long long entries = chain->GetEntries();
// 	long long entry100 = entries / 100;
// 	printf("Filling front back energy graph   0%%");
//     fflush(stdout);

// 	// points on the correlated energy graph
// 	int points = 0;
// 	// loop entries and generate correlated energy graph
// 	for (long long  entry = 0; entry < entries; ++entry)
// 	{
// 		// show process
// 		if (entry % entry100 == 0)
// 		{
//             printf("\b\b\b\b%3lld%%", entry / entry100);
//             fflush(stdout);
//         }

// 		// read event
// 		chain->GetEntry(entry);

// 		// ignore the events with hit over 1 now
// 		if (correlation_.front_hit != 1 || correlation_.back_hit != 1)
// 		{
// 			continue;
// 		}

// 		// use single strip as reference
// 		if (correlation_.back_strip[0] != ref_back) continue;

// 		// check wheter the correlated event energy is suitable for normalize
// 		if (!NormalizeFrontEnergyCheck(correlation_, iteration)) continue;

// 		// add point to the graph if appropriate
// 		g_front_back_energy[correlation_.front_strip[0]]->SetPoint
// 		(
// 			points++,
// 			correlation_.front_energy[0],
// 			correlation_.back_energy[0]
// 		);
// 	}
// 	// show finish
// 	printf("\b\b\b\b100%%\n");

// 	// fitting now
// 	std::cout << "front strip parameters\n";
// 	for (size_t i = 0; i < FrontStrip(); ++i)
// 	{
// 		// only fits when over 10 points
// 		if (g_front_back_energy[i]->GetN() > 10)
// 		{
// 			// fitting function
// 			TF1 *energy_fit = new TF1("efit", "pol1", 0, 60000);
// 			// set initial value
// 			energy_fit->SetParameter(0, 0.0);
// 			energy_fit->SetParameter(1, 1.0);
// 			// fit
// 			g_front_back_energy[i]->Fit(energy_fit, "QR+");
// 			// store the normalized parameters
// 			normalize_front_param_[i][0] = energy_fit->GetParameter(0);
// 			normalize_front_param_[i][1] = energy_fit->GetParameter(1);
// 		}
// 		// store the graph
// 		g_front_back_energy[i]->Write(TString::Format("gf%ld", i));
// 		// print normalized paramters on screen
// 		std::cout << i << " " << normalize_front_param_[i][0]
// 			<< ", " << normalize_front_param_[i][1] << "\n";
// 	}

// 	return;
// }


// void Detector::NormalizeSecondSide(
// 	TTree *chain,
// 	unsigned short ref_front,
// 	unsigned short,
// 	TGraph **,
// 	TGraph **g_back_front_energy,
// 	bool iteration
// ) {
// 	// show process
// 	long long entries = chain->GetEntries();
// 	long long entry100 = entries / 100;
// 	printf("Filling back front energy graph   0%%");
//     fflush(stdout);

// 	// points on correlated energy graph
// 	int points = 0;
// 	// loop entries and generate correlated energy graph
// 	for (long long entry = 0; entry < entries; ++entry)
// 	{
// 		// show process
// 		if (entry % entry100 == 0)
// 		{
//             printf("\b\b\b\b%3lld%%", entry / entry100);
//             fflush(stdout);
//         }

// 		// read event
// 		chain->GetEntry(entry);

// 		// hit over 1, ignore it now
// 		if (correlation_.front_hit != 1 || correlation_.back_hit != 1)
// 			continue;

// 		// use single strip at the front side
// 		if (correlation_.front_strip[0] != ref_front) continue;

// 		// check wheter the correlated event energy is suitable for normalize
// 		if (!NormalizeBackEnergyCheck(correlation_, iteration)) continue;

// 		// fill the appropiate point into the graph
// 		g_back_front_energy[correlation_.back_strip[0]]->SetPoint(
// 			points++,
// 			correlation_.back_energy[0],
// 			NormalizeEnergy(
// 				0,
// 				correlation_.front_strip[0],
// 				correlation_.front_energy[0]
// 			)
// 		);
// 	}
// 	// show process
// 	printf("\b\b\b\b100%%\n");

// 	// fitting now
// 	std::cout << "back strip parameters\n";
// 	for (size_t i = 0; i < BackStrip(); ++i) {
// 		// fit only if over 5 points
// 		if (g_back_front_energy[i]->GetN() > 10) {
// 			// fitting function
// 			TF1 *energy_fit = new TF1("efit", "pol1", 0, 60000);
// 			// set initial value
// 			energy_fit->SetParameter(0, 0.0);
// 			energy_fit->SetParameter(1, 1.0);
// 			// fit
// 			g_back_front_energy[i]->Fit(energy_fit, "QR+");
// 			// store normalized paramters
// 			normalize_back_param_[i][0] = energy_fit->GetParameter(0);
// 			normalize_back_param_[i][1] = energy_fit->GetParameter(1);
// 		}
// 		// store the graph
// 		g_back_front_energy[i]->Write(TString::Format("gb%ld", i));
// 		// print the normalized parameters
// 		std::cout << i << " " << normalize_back_param_[i][0]
// 			<< ", " << normalize_back_param_[i][1] << "\n";
// 	}
// 	return;
// }


// int Detector::Normalize
// (
// 	int length,
// 	unsigned short ref_front,
// 	unsigned short ref_back,
// 	bool iteration
// )
// {
// 	// setup input chain
// 	TChain *chain = new TChain("tree", "chain");
// 	for (int i = 0; i < length; ++i) {
// 		TString file_name;
// 		file_name.Form(
// 			"%s%s%s-correlated-%04d.root",
// 			kGenerateDataPath, kCorrelationDir, name_.c_str(), run_
// 		);
// 		chain->AddFile(file_name);
// 	}
// 	// setup input correlated tree
// 	SetupInputCorrelatedTree(chain);

// 	// initialize normalized parameters
// 	for (size_t i = 0; i < FrontStrip(); ++i) {
// 		normalize_front_param_[i][0] = 0.0;
// 		normalize_front_param_[i][1] = 1.0;
// 	}
// 	for (size_t i = 0; i < BackStrip(); ++i) {
// 		normalize_back_param_[i][0] = 0.0;
// 		normalize_back_param_[i][1] = 1.0;
// 	}

// 	// read normalized parameters from file if in iteration mode
// 	if (iteration && ReadNormalizeParameters()) {
// 		std::cerr << "Error: read normalize parameters from file failed.\n";
// 		return -1;
// 	}

// 	// setup normalize record root file
// 	TString normalize_file_name;
// 	normalize_file_name.Form
// 	(
// 		"%s%s%s-normalize-%04d.root",
// 		kGenerateDataPath, kNormalizeDir, name_.c_str(), run_
// 	);
// 	TFile *opf = new TFile(normalize_file_name, "recreate");

// 	TGraph *g_front_back_energy[64];
// 	TGraph *g_back_front_energy[64];
// 	for (size_t i = 0; i < FrontStrip(); ++i)
// 	{
// 		g_front_back_energy[i] = new TGraph;
// 	}
// 	for (size_t i = 0; i < BackStrip(); ++i)
// 	{
// 		g_back_front_energy[i] = new TGraph;
// 	}

// 	// normalize the first side
// 	NormalizeFirstSide(
// 		chain,
// 		ref_front,
// 		ref_back,
// 		g_front_back_energy,
// 		g_back_front_energy,
// 		iteration
// 	);
// 	// normalize the second side
// 	NormalizeSecondSide(
// 		chain,
// 		ref_front,
// 		ref_back,
// 		g_front_back_energy,
// 		g_back_front_energy,
// 		iteration
// 	);

// 	opf->Close();

// 	if (WriteNormalizeParameters())
// 	{
// 		std::cerr << "Error: write normalize paramters to file failed.\n";
// 		return -1;
// 	}

// 	return 0;
// }




// int Detector::NormalCorrelate() {
// 	// setup input file
// 	TString input_file_name;
// 	input_file_name.Form(
// 		"%s%s%s-correlated-%04d.root",
// 		kGenerateDataPath, kCorrelationDir, name_.c_str(), run_
// 	);
// 	TFile *ipf = new TFile(input_file_name, "read");
// 	TTree *ipt = (TTree*)ipf->Get("tree");
// 	if (!ipt) {
// 		std::cerr << "Error: get tree from "
// 			<< input_file_name << " failed.\n";
// 		return -1;
// 	}
// 	// setup input tree
// 	SetupInputCorrelatedTree(ipt);

// 	// setup output file
// 	TString output_file_name;
// 	output_file_name.Form(
// 		"%s%s%s-normalized-correlated-%04d.root",
// 		kGenerateDataPath, kCorrelationDir, name_.c_str(), run_
// 	);
// 	TFile *opf = new TFile(output_file_name, "recreate");
// 	TTree *opt = new TTree("tree", "tree after normalize and correlation");
// 	// front hit
// 	unsigned short front_hit;
// 	// back hit
// 	unsigned short back_hit;
// 	// normalized front energy
// 	double front_norm_energy[8];
// 	// normalized back energy
// 	double back_norm_energy[8];
// 	// setup output tree
// 	opt->Branch("front_hit", &front_hit, "fhit/s");
// 	opt->Branch("front_energy", front_norm_energy, "fe[fhit]/D");
// 	opt->Branch("back_hit", &back_hit, "bhit/s");
// 	opt->Branch("back_energy", back_norm_energy, "be[bhit]/D");

// 	// read normalize paramters
// 	if (ReadNormalizeParameters()) {
// 		return -1;
// 	}

// 	// show process
// 	printf("building normalized result   0%%");
// 	fflush(stdout);
// 	long long entries = ipt->GetEntries();
// 	long long entry100 = entries / 100;
// 	for (long long entry = 0; entry < entries; ++entry) {
// 		// show process
// 		if (entry % entry100 == 0) {
// 			printf("\b\b\b\b%3lld%%", entry / entry100);
// 			fflush(stdout);
// 		}

// 		ipt->GetEntry(entry);
// 		front_hit = correlation_.front_hit;
// 		back_hit = correlation_.back_hit;
// 		if (front_hit < 8 && back_hit < 8) {
// 			for (unsigned short i = 0; i < front_hit; ++i) {
// 				front_norm_energy[i] =
// 					NormalizeEnergy(0, correlation_.front_strip[i], correlation_.front_energy[i]);
// 			}
// 			for (unsigned short i = 0; i < back_hit; ++i) {
// 				back_norm_energy[i] =
// 					NormalizeEnergy(1, correlation_.back_strip[i], correlation_.back_energy[i]);
// 			}
// 		}
// 		opt->Fill();
// 	}
// 	// show finish
// 	printf("\b\b\b\b100%%\n");

// 	opt->Write();
// 	opf->Close();

// 	ipf->Close();

// 	return 0;
// }



// //-----------------------------------------------------------------------------
// //										merge
// //-----------------------------------------------------------------------------


// void Detector::SetupMergedOutput(TFile **output_file, TTree **output_tree)
// {
// 	// merged event file name
// 	TString merged_file_name;
// 	merged_file_name.Form(
// 		"%s%s%s-merged-%04d.root",
// 		kGenerateDataPath, kMergedDir, name_.c_str(), run_
// 	);
// 	// setup merged file
// 	TFile *merged_file = new TFile(merged_file_name, "recreate");
// 	// setup merged tree
// 	TTree *merged_tree = new TTree(
// 		"tree", TString::Format("merged events of %s", name_.c_str())
// 	);
// 	// setup merged tree branches
// 	merged_tree->Branch("timestamp", &merged_.timestamp, "ts/L");
// 	merged_tree->Branch("hit", &merged_.hit, "hit/s");
// 	merged_tree->Branch("front_strip", merged_.front_strip, "fs[hit]/D");
// 	merged_tree->Branch("back_strip", merged_.back_strip, "bs[hit]/D");
// 	merged_tree->Branch("time", merged_.time, "t[hit]/D");
// 	merged_tree->Branch("energy", merged_.energy, "e[hit]/D");

// 	// setup output
// 	*output_file = merged_file;
// 	*output_tree = merged_tree;

// 	return;
// }


// void Detector::SetupUnmergedOutput
// (
// 	TFile **output_file,
// 	TTree **output_tree
// )
// {
// 	// setup residual correlated file
// 	TString file_name;
// 	file_name.Form(
// 		"%s%s%s-unmerged-%04d.root",
// 		kGenerateDataPath, kMergedDir, name_.c_str(), run_
// 	);
// 	*output_file = new TFile(file_name, "recreate");
// 	// setup output tree
// 	*output_tree = new TTree
// 	(
// 		"tree",
// 		TString::Format("unmerged tree of %s", name_.c_str())
// 	);
// 	// setup output data
// 	SetupOutputCorrelatedTree(*output_tree);

// 	return;
// }



// inline bool EnergyCut
// (
// 	double front_energy,
// 	double back_energy,
// 	double energy_cut
// )
// {
// 	return
// 		std::abs((front_energy - back_energy) / (front_energy + back_energy))
// 			< energy_cut;
// }


// int Detector::Merge(double energy_cut)
// {
// 	// setup input file
// 	TString correlated_file_name;
// 	correlated_file_name.Form(
// 		"%s%s%s-correlated-%04d.root",
// 		kGenerateDataPath, kCorrelationDir, name_.c_str(), run_
// 	);
// 	TFile *correlated_file = new TFile(correlated_file_name, "read");
// 	TTree *correlated_tree = (TTree*)correlated_file->Get("tree");
// 	if (!correlated_tree)
// 	{
// 		std::cerr << "Error: get tree from " << correlated_file_name << " failed\n";
// 		return -1;
// 	}
// 	SetupInputCorrelatedTree(correlated_tree);

// 	// setup output file for merged events
// 	TFile *merged_file = nullptr;
// 	TTree *merged_tree = nullptr;
// 	SetupMergedOutput(&merged_file, &merged_tree);

// 	// setup output file for unmerged events
// 	TFile *unmerged_file = nullptr;
// 	TTree *unmerged_tree = nullptr;
// 	SetupUnmergedOutput(&unmerged_file, &unmerged_tree);

// 	// read normalized parameters from file
// 	if (ReadNormalizeParameters())
// 	{
// 		std::cerr << "Error: Read normalize parameters failed.\n";
// 		return -1;
// 	}

// 	// show process of merging events
// 	printf("Merging %s   0%%", name_.c_str());
// 	fflush(stdout);
// 	// entry number
// 	long long entries = correlated_tree->GetEntries();
// 	// 1/100 of entry
// 	long long entry100 = entries / 100;
// 	// merged events number
// 	long long merged_events = 0;

// 	// loop for merging events
// 	for (long long entry = 0; entry < entries; ++entry)
// 	{
// 		// show process
// 		if (entry % entry100 == 0)
// 		{
// 			printf("\b\b\b\b%3lld%%", entry / entry100);
// 			fflush(stdout);
// 		}
// 		correlated_tree->GetEntry(entry);

// 		merged_.timestamp = correlation_.timestamp;
// 		merged_.hit = 0;
// 		double front_total_energy = 0.0;
// 		double back_total_energy = 0.0;

// 		// calculate normalied energy and total energy for front side
// 		for (unsigned short i = 0; i < correlation_.front_hit; ++i)
// 		{
// 			correlation_.front_energy[i] = NormalizeEnergy
// 				(
// 					0,
// 					correlation_.front_strip[i],
// 					correlation_.front_energy[i]
// 				);
// 			front_total_energy += correlation_.front_energy[i];
// 		}
// 		// calculate normalized energy and total energy for back side
// 		for (unsigned short i = 0; i < correlation_.back_hit; ++i)
// 		{
// 			correlation_.back_energy[i] = NormalizeEnergy
// 				(
// 					1,
// 					correlation_.back_strip[i],
// 					correlation_.back_energy[i]
// 				);
// 			back_total_energy += correlation_.back_energy[i];
// 		}

// 		if (correlation_.back_hit == 1)
// 		{

// 			if (correlation_.front_hit == 1)
// 			{

// 				if (EnergyCut(front_total_energy, back_total_energy, energy_cut))
// 				{
// 					++merged_events;
// 					merged_.hit = 1;
// 					merged_.front_strip[0] = correlation_.front_strip[0];
// 					merged_.back_strip[0] = correlation_.front_strip[0];
// 					merged_.time[0] = correlation_.back_time[0];
// 					merged_.energy[0] = back_total_energy;
// 				}

// 			}
// 			else if (correlation_.front_hit == 2)
// 			{

// 				if
// 				(
// 					EnergyCut(front_total_energy, back_total_energy, energy_cut)
// 					&& correlation_.front_strip[0] + 1 == correlation_.front_strip[1]
// 				)
// 				{

// 					++merged_events;
// 					merged_.hit = 1;
// 					merged_.front_strip[0] = correlation_.front_strip[0]
// 						+ correlation_.front_energy[1] / front_total_energy;
// 					merged_.back_strip[0] = correlation_.back_strip[0];
// 					merged_.time[0] = correlation_.back_time[0];
// 					merged_.energy[0] = back_total_energy;

// 				}

// 			}
// 			else if(correlation_.front_hit == 3)
// 			{

// 				if
// 				(
// 					EnergyCut(front_total_energy, back_total_energy, energy_cut)
// 					&& correlation_.front_strip[0] + correlation_.front_strip[2]
// 						== correlation_.front_strip[1] * 2
// 				)
// 				{

// 					++merged_events;
// 					merged_.hit = 1;
// 					merged_.front_strip[0] = correlation_.front_strip[1];
// 					merged_.back_strip[0] = correlation_.back_strip[0];
// 					merged_.time[0] = correlation_.back_time[0];
// 					merged_.energy[0] = back_total_energy;

// 				}

// 			}

// 		}
// 		else if (correlation_.back_hit == 2)
// 		{

// 			if (correlation_.front_hit == 1)
// 			{

// 				if (
// 					EnergyCut(front_total_energy, back_total_energy, energy_cut)
// 					&& correlation_.back_strip[0] + 1 == correlation_.back_strip[1]
// 				)
// 				{

// 					++merged_events;
// 					merged_.hit = 1;
// 					merged_.front_strip[0] = correlation_.front_strip[0];
// 					merged_.back_strip[0] = correlation_.back_strip[0]
// 						+ correlation_.back_energy[1] / back_total_energy;
// 					merged_.time[0] = correlation_.front_time[0];
// 					merged_.energy[0] = front_total_energy;

// 				}

// 			} else if (correlation_.front_hit == 2) {

// 				if
// 				(
// 					EnergyCut(
// 						correlation_.front_energy[0],
// 						correlation_.back_energy[0],
// 						energy_cut
// 					)
// 					&& EnergyCut(
// 						correlation_.front_energy[1],
// 						correlation_.back_energy[1],
// 						energy_cut
// 					)
// 				)
// 				{
// 					++merged_events;
// 					merged_.hit = 2;
// 					if (correlation_.front_energy[0] > correlation_.front_energy[1])
// 					{
// 						merged_.front_strip[0] = correlation_.front_strip[0];
// 						merged_.front_strip[1] = correlation_.front_strip[1];
// 						merged_.back_strip[0] = correlation_.back_strip[0];
// 						merged_.back_strip[1] = correlation_.back_strip[1];
// 						merged_.time[0] = correlation_.front_time[0];
// 						merged_.time[1] = correlation_.front_time[1];
// 						merged_.energy[0] = correlation_.front_energy[0];
// 						merged_.energy[1] = correlation_.front_energy[1];
// 					}
// 					else
// 					{
// 						merged_.front_strip[0] = correlation_.front_strip[1];
// 						merged_.front_strip[1] = correlation_.front_strip[0];
// 						merged_.back_strip[0] = correlation_.back_strip[1];
// 						merged_.back_strip[1] = correlation_.back_strip[0];
// 						merged_.time[0] = correlation_.front_time[1];
// 						merged_.time[1] = correlation_.front_time[0];
// 						merged_.energy[0] = correlation_.front_energy[1];
// 						merged_.energy[1] = correlation_.front_energy[0];
// 					}

// 				}
// 				else if
// 				(
// 					EnergyCut(
// 						correlation_.front_energy[0],
// 						correlation_.back_energy[1],
// 						energy_cut
// 					)
// 					&& EnergyCut(
// 						correlation_.front_energy[1],
// 						correlation_.back_energy[0],
// 						energy_cut
// 					)
// 				)
// 				{

// 					++merged_events;
// 					merged_.hit = 2;
// 					if (correlation_.front_energy[0] > correlation_.front_energy[1])
// 					{
// 						merged_.front_strip[0] = correlation_.front_strip[0];
// 						merged_.front_strip[1] = correlation_.front_strip[1];
// 						merged_.back_strip[0] = correlation_.back_strip[1];
// 						merged_.back_strip[1] = correlation_.back_strip[0];
// 						merged_.time[0] = correlation_.front_time[0];
// 						merged_.time[1] = correlation_.front_time[1];
// 						merged_.energy[0] = correlation_.front_energy[0];
// 						merged_.energy[1] = correlation_.front_energy[1];
// 					}
// 					else
// 					{
// 						merged_.front_strip[0] = correlation_.front_strip[1];
// 						merged_.front_strip[1] = correlation_.front_strip[0];
// 						merged_.back_strip[0] = correlation_.back_strip[0];
// 						merged_.back_strip[1] = correlation_.back_strip[1];
// 						merged_.time[0] = correlation_.front_time[1];
// 						merged_.time[1] = correlation_.front_time[0];
// 						merged_.energy[0] = correlation_.front_energy[1];
// 						merged_.energy[1] = correlation_.front_energy[0];
// 					}

// 				}
// 				else if
// 				(
// 					EnergyCut(front_total_energy, back_total_energy, energy_cut)
// 					&& correlation_.front_strip[0] + 1 == correlation_.front_strip[1]
// 					&& correlation_.back_strip[0] + 1 == correlation_.back_strip[1]
// 				)
// 				{
// 					++merged_events;
// 					merged_.hit = 1;
// 					merged_.front_strip[0] = correlation_.front_strip[0]
// 						+ correlation_.front_energy[1] / front_total_energy;
// 					merged_.back_strip[0] = correlation_.back_strip[0]
// 						+ correlation_.back_energy[1] / back_total_energy;
// 					merged_.time[0] = correlation_.back_time[0];
// 					merged_.energy[0] = back_total_energy;
// 				}

// 			}
// 			else if (correlation_.front_hit == 3)
// 			{

// 				if
// 				(
// 					EnergyCut
// 					(
// 						correlation_.back_energy[0],
// 						correlation_.front_energy[0] + correlation_.front_energy[1],
// 						energy_cut
// 					)
// 					&& EnergyCut
// 					(
// 						correlation_.back_energy[1],
// 						correlation_.front_energy[2],
// 						energy_cut
// 					)
// 					&& correlation_.front_strip[0] + 1 == correlation_.front_strip[1]
// 				)
// 				{

// 					++merged_events;
// 					merged_.hit = 2;
// 					if (correlation_.back_energy[0] > correlation_.back_energy[1])
// 					{
// 						merged_.front_strip[0] = correlation_.front_strip[0]
// 							+ correlation_.front_energy[1] / front_total_energy;
// 						merged_.front_strip[1] = correlation_.front_strip[2];
// 						merged_.back_strip[0] = correlation_.back_strip[0];
// 						merged_.back_strip[1] = correlation_.back_strip[1];
// 						merged_.time[0] = correlation_.back_time[0];
// 						merged_.time[1] = correlation_.back_time[1];
// 						merged_.energy[0] = correlation_.back_energy[0];
// 						merged_.energy[1] = correlation_.back_energy[1];
// 					}
// 					else
// 					{
// 						merged_.front_strip[0] = correlation_.front_strip[2];
// 						merged_.front_strip[1] = correlation_.front_strip[0]
// 							+ correlation_.front_energy[1] / front_total_energy;
// 						merged_.back_strip[0] = correlation_.back_strip[1];
// 						merged_.back_strip[1] = correlation_.back_strip[0];
// 						merged_.time[0] = correlation_.back_time[1];
// 						merged_.time[1] = correlation_.back_time[0];
// 						merged_.energy[0] = correlation_.back_energy[1];
// 						merged_.energy[1] = correlation_.back_energy[0];
// 					}

// 				}
// 				else if
// 				(
// 					EnergyCut
// 					(
// 						correlation_.back_energy[0],
// 						correlation_.front_energy[2],
// 						energy_cut
// 					)
// 					&& EnergyCut
// 					(
// 						correlation_.back_energy[1],
// 						correlation_.front_energy[0] + correlation_.front_energy[1],
// 						energy_cut
// 					)
// 					&& correlation_.front_strip[0] + 1 == correlation_.front_strip[1]
// 				)
// 				{

// 					++merged_events;
// 					merged_.hit = 2;
// 					if (correlation_.back_energy[0] > correlation_.back_energy[1])
// 					{
// 						merged_.front_strip[0] = correlation_.front_strip[2];
// 						merged_.front_strip[1] = correlation_.front_strip[0]
// 							+ correlation_.front_energy[1] / front_total_energy;
// 						merged_.back_strip[0] = correlation_.back_strip[0];
// 						merged_.back_strip[1] = correlation_.back_strip[1];
// 						merged_.time[0] = correlation_.back_time[0];
// 						merged_.time[1] = correlation_.back_time[1];
// 						merged_.energy[0] = correlation_.back_energy[0];
// 						merged_.energy[1] = correlation_.back_energy[1];
// 					}
// 					else
// 					{
// 						merged_.front_strip[0] = correlation_.front_strip[0]
// 							+ correlation_.front_energy[1] / front_total_energy;
// 						merged_.front_strip[1] = correlation_.front_strip[2];
// 						merged_.back_strip[0] = correlation_.back_strip[1];
// 						merged_.back_strip[1] = correlation_.back_strip[0];
// 						merged_.time[0] = correlation_.back_time[1];
// 						merged_.time[1] = correlation_.back_time[0];
// 						merged_.energy[0] = correlation_.back_energy[1];
// 						merged_.energy[1] = correlation_.back_energy[0];
// 					}

// 				}
// 				else if
// 				(
// 					EnergyCut
// 					(
// 						correlation_.back_energy[0],
// 						correlation_.front_energy[1] + correlation_.front_energy[2],
// 						energy_cut
// 					)
// 					&& EnergyCut
// 					(
// 						correlation_.back_energy[1],
// 						correlation_.front_energy[0],
// 						energy_cut
// 					)
// 					&& correlation_.front_strip[1] + 1 == correlation_.front_strip[2]
// 				)
// 				{

// 					++merged_events;
// 					merged_.hit = 2;
// 					if (correlation_.back_energy[0] > correlation_.back_energy[1])
// 					{
// 						merged_.front_strip[0] = correlation_.front_strip[1]
// 							+ correlation_.front_energy[2] / front_total_energy;
// 						merged_.front_strip[1] = correlation_.front_strip[0];
// 						merged_.back_strip[0] = correlation_.back_strip[0];
// 						merged_.back_strip[1] = correlation_.back_strip[1];
// 						merged_.time[0] = correlation_.back_time[0];
// 						merged_.time[1] = correlation_.back_time[1];
// 						merged_.energy[0] = correlation_.back_energy[0];
// 						merged_.energy[1] = correlation_.back_energy[1];
// 					}
// 					else
// 					{
// 						merged_.front_strip[0] = correlation_.front_strip[0];
// 						merged_.front_strip[1] = correlation_.front_strip[1]
// 							+ correlation_.front_energy[2] / front_total_energy;
// 						merged_.back_strip[0] = correlation_.back_strip[1];
// 						merged_.back_strip[1] = correlation_.back_strip[0];
// 						merged_.time[0] = correlation_.back_time[1];
// 						merged_.time[1] = correlation_.back_time[0];
// 						merged_.energy[0] = correlation_.back_energy[1];
// 						merged_.energy[1] = correlation_.back_energy[0];
// 					}

// 				}
// 				else if
// 				(
// 					EnergyCut
// 					(
// 						correlation_.back_energy[0],
// 						correlation_.front_energy[0],
// 						energy_cut
// 					)
// 					&& EnergyCut
// 					(
// 						correlation_.back_energy[1],
// 						correlation_.front_energy[1] + correlation_.front_energy[2],
// 						energy_cut
// 					)
// 					&& correlation_.front_strip[1] + 1 == correlation_.front_strip[2]
// 				)
// 				{

// 					++merged_events;
// 					merged_.hit = 2;
// 					if (correlation_.back_energy[0] > correlation_.back_energy[1])
// 					{
// 						merged_.front_strip[0] = correlation_.front_strip[0];
// 						merged_.front_strip[1] = correlation_.front_strip[1]
// 							+ correlation_.front_energy[2] / front_total_energy;
// 						merged_.back_strip[0] = correlation_.back_strip[0];
// 						merged_.back_strip[1] = correlation_.back_strip[1];
// 						merged_.time[0] = correlation_.back_time[0];
// 						merged_.time[1] = correlation_.back_time[1];
// 						merged_.energy[0] = correlation_.back_energy[0];
// 						merged_.energy[1] = correlation_.back_energy[1];
// 					}
// 					else
// 					{
// 						merged_.front_strip[0] = correlation_.front_strip[1]
// 							+ correlation_.front_energy[2] / front_total_energy;
// 						merged_.front_strip[1] = correlation_.front_strip[0];
// 						merged_.back_strip[0] = correlation_.back_strip[1];
// 						merged_.back_strip[1] = correlation_.back_strip[0];
// 						merged_.time[0] = correlation_.back_time[1];
// 						merged_.time[1] = correlation_.back_time[0];
// 						merged_.energy[0] = correlation_.back_energy[1];
// 						merged_.energy[1] = correlation_.back_energy[0];
// 					}

// 				// } else if (
// 				// 	EnergyCut(front_total_energy, back_total_energy, energy_cut)
// 				// 	&& correlation_.back_strip[0] + 1 == correlation_.back_strip[1]
// 				// 	&& correlation_.front_strip[0] + correlation_.front_strip[2]
// 				// 		== correlation_.back_strip[1] * 2
// 				// ) {

// 				// 	++merged_events;
// 				// 	merged_.hit = 1;
// 				// 	merged_.front_strip[0] = correlation_.front_strip[1];
// 				// 	merged_.back_strip[0] = correlation_.back_strip[0] +
// 				// 		correlation_.back_energy[1] / back_total_energy;
// 				// 	merged_.time[0] = correlation_.back_time[0];
// 				//  merged_.energy[0] = back_total_energy;

// 				}


// 			}
// 			else if (correlation_.front_hit == 4)
// 			{
// 			}

// 		}

// 		merged_file->cd();
// 		merged_tree->Fill();

// 		if (!merged_.hit) {
// 			unmerged_file->cd();
// 			unmerged_tree->Fill();
// 		}
// 	}
// 	// show finished
// 	printf("\b\b\b\b100%%\n");

// 	// write and close output file
// 	merged_file->cd();
// 	merged_tree->Write();
// 	merged_file->Close();

// 	// write residual correlated output file
// 	unmerged_file->cd();
// 	unmerged_tree->Write();
// 	unmerged_file->Close();

// 	// close input file
// 	correlated_file->Close();

// 	std::cout << "merged rate " << merged_events << " / " << entries
// 		<< " " << double(merged_events) / entries << "\n";

// 	return 0;
// }

// //-----------------------------------------------------------------------------
// // 									telescope
// //-----------------------------------------------------------------------------

// int Detector::ReadMergedEvents()
// {
// 	// setup input file
// 	TString file_name;
// 	file_name.Form(
// 		"%s%s%s-merged-%04d.root",
// 		kGenerateDataPath, kMergedDir, name_.c_str(), run_
// 	);
// 	input_merged_file_ = new TFile(file_name, "read");
// 	// setup input tree
// 	input_merged_tree_ = (TTree*)input_merged_file_->Get("tree");
// 	// check input tree
// 	if (!input_merged_tree_)
// 	{
// 		std::cerr << "Error: get tree from " << file_name << " failed.\n";
// 		return -1;
// 	}
// 	// setup branches
// 	input_merged_tree_->SetBranchAddress("hit", &merged_.hit);
// 	input_merged_tree_->SetBranchAddress("timestamp", &merged_.timestamp);
// 	input_merged_tree_->SetBranchAddress("front_strip", merged_.front_strip);
// 	input_merged_tree_->SetBranchAddress("back_strip", merged_.back_strip);
// 	input_merged_tree_->SetBranchAddress("time", merged_.time);
// 	input_merged_tree_->SetBranchAddress("energy", merged_.energy);

// 	merged_events_ = input_merged_tree_->GetEntries();
// 	return 0;
// }


// std::shared_ptr<Detector> CreateDetector(int run, const std::string &name) {
// 	if (name == "t0d1")
// 	{
// 		return std::make_shared<T0D1>(run);
// 	}
// 	else if (name == "t0d2")
// 	{
// 		return std::make_shared<T0D2>(run);
// 	}
// 	else if (name == "t0d3")
// 	{
// 		return std::make_shared<T0D3>(run);
// 	}
// 	return nullptr;
// }



// //-----------------------------------------------------------------------------
// // 									DSSD
// //-----------------------------------------------------------------------------

// DSSD::DSSD(int run, const std::string &name)
// : Detector(run, name) {
// }


// void DSSD::SetupOutputCorrelatedTree(TTree *tree) {
// 	tree->Branch("timestamp", &correlation_.timestamp, "ts/L");
// 	tree->Branch("front_hit", &correlation_.front_hit, "fhit/s");
// 	tree->Branch("back_hit", &correlation_.back_hit, "bhit/s");
// 	tree->Branch("front_strip", correlation_.front_strip, "fs[fhit]/s");
// 	tree->Branch("back_strip", correlation_.back_strip, "bs[bhit]/s");
// 	tree->Branch("front_time", correlation_.front_time, "ft[fhit]/D");
// 	tree->Branch("back_time", correlation_.back_time, "bt[bhit]/D");
// 	tree->Branch("front_energy", correlation_.front_energy, "fe[fhit]/D");
// 	tree->Branch("back_energy", correlation_.back_energy, "be[bhit]/D");
// }


// void DSSD::SetupInputCorrelatedTree(TTree *tree) {
// 	tree->SetBranchAddress("timestamp", &correlation_.timestamp);
// 	tree->SetBranchAddress("front_hit", &correlation_.front_hit);
// 	tree->SetBranchAddress("back_hit", &correlation_.back_hit);
// 	tree->SetBranchAddress("front_strip", correlation_.front_strip);
// 	tree->SetBranchAddress("back_strip", correlation_.back_strip);
// 	tree->SetBranchAddress("front_time", correlation_.front_time);
// 	tree->SetBranchAddress("back_time", correlation_.back_time);
// 	tree->SetBranchAddress("front_energy", correlation_.front_energy);
// 	tree->SetBranchAddress("back_energy", correlation_.back_energy);
// }


// //-----------------------------------------------------------------------------
// // 									T0D1
// //-----------------------------------------------------------------------------

// T0D1::T0D1(int run)
// : DSSD(run, "t0d1") {
// }


// bool T0D1::NormalizeFrontEnergyCheck(
// 	const CorrelatedEvent &correlation,
// 	bool iteration
// )
// {
// 	if
// 	(
// 		correlation.front_energy[0] < 5000
// 		|| correlation.back_energy[0] < 5000
// 	)
// 		return false;


// 	if (!iteration)
// 	{
// 		if (correlation.front_strip[0] < 16)
// 		{
// 		}
// 		else if (correlation.front_strip[0] <= 32)
// 		{
// 			if
// 			(
// 				(
// 					correlation.front_energy[0] > 18000
// 					&& correlation.front_energy[0] < 30000
// 				)
// 				||
// 				(
// 					correlation.back_energy[0] > 30000
// 					&& correlation.back_energy[0] < 44000
// 				)
// 			)
// 			{
// 				return false;
// 			}

// 		}
// 		else if (correlation.front_strip[0] < 48)
// 		{
// 			if
// 			(
// 				(
// 					correlation.front_energy[0] > 23000
// 					&& correlation.front_energy[0] < 27000
// 				)
// 				||
// 				(
// 					correlation.back_energy[0] > 34000
// 					&& correlation.back_energy[0] < 41000
// 				)
// 			)
// 			{
// 				return false;
// 			}

// 			// if (
// 			// 	abs(correlation.front_energy[0] - norm_back_energy)
// 			// 	> 0.6 * norm_back_energy
// 			// ) return false;

// 		}
// 		else
// 		{

// 		}

// 	}
// 	else
// 	{
// 		if
// 		(
// 			abs
// 			(
// 				NormalizeEnergy
// 				(
// 					0,
// 					correlation.front_strip[0],
// 					correlation.front_energy[0]
// 				)
// 				- correlation.back_energy[0]
// 			)
// 			> 1000
// 		) return false;

// 		if (correlation.front_strip[0] < 16)
// 		{


// 		}
// 		else if (correlation.front_strip[0] <= 32)
// 		{
// 			if
// 			(
// 				(
// 					correlation.front_energy[0] > 23000
// 					&& correlation.front_energy[0] < 26000
// 				)
// 				||
// 				(
// 					correlation.back_energy[0] > 34000
// 					&& correlation.back_energy[0] < 40000
// 				)
// 			)
// 			{
// 				return false;
// 			}


// 		}
// 		else if (correlation.front_strip[0] < 48)
// 		{
// 			if
// 			(
// 				(
// 					correlation.front_energy[0] > 17000
// 					&& correlation.front_energy[0] < 35000
// 				)
// 				||
// 				(
// 					correlation.back_energy[0] > 25000
// 					&& correlation.back_energy[0] < 49000
// 				)
// 			)
// 			{
// 				return false;
// 			}

// 		}
// 		else
// 		{
// 		}
// 	}

// 	return true;
// }



// bool T0D1::NormalizeBackEnergyCheck(const CorrelatedEvent &correlation, bool iteration) {
// 	if
// 	(
// 		correlation.front_energy[0] < 5000
// 		|| correlation.back_energy[0] < 5000
// 	)
// 		return false;

// 	if (!iteration)
// 	{
// 		// if (correlation.back_strip[0] < 16) {
// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 20000 && correlation.back_energy[0] < 22000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 30000 && correlation.back_energy[0] > 34000
// 		// 		)
// 		// 	)) return false;
// 		// } else if (correlation.back_strip[0] < 32) {
// 		// 	if (
// 		// 		abs(correlation.front_energy[0] - correlation.back_energy[0])
// 		// 		> 0.3 * correlation.back_energy[0]
// 		// 	) return false;

// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 20000 && correlation.back_energy[0] < 24000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 30000 && correlation.back_energy[0] > 38000
// 		// 		)
// 		// 	)) return false;
// 		// } else if (correlation.back_strip[0] < 48) {
// 		// 	if (
// 		// 		abs(correlation.front_energy[0] - correlation.back_energy[0])
// 		// 		> 0.5 * correlation.back_energy[0]
// 		// 	) return false;

// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 19000 && correlation.back_energy[0] < 28000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 29000 && correlation.back_energy[0] > 42000
// 		// 		)
// 		// 	)) return false;

// 		// } else {

// 		// }

// 	}
// 	else
// 	{
// 		if
// 		(
// 			abs(
// 				NormalizeEnergy(
// 					0,
// 					correlation.front_strip[0],
// 					correlation.front_energy[0]
// 				)
// 				- NormalizeEnergy(
// 					1,
// 					correlation.back_strip[0],
// 					correlation.back_energy[0]
// 				)
// 			)
// 			> 1000
// 		)
// 			return false;


// 		// if (correlation.back_strip[0] < 16) {
// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 20000 && correlation.back_energy[0] < 22000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 30000 && correlation.back_energy[0] > 34000
// 		// 		)
// 		// 	)) return false;
// 		// } else if (correlation.back_strip[0] < 32) {
// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 20000 && correlation.back_energy[0] < 24000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 30000 && correlation.back_energy[0] > 38000
// 		// 		)
// 		// 	)) return false;
// 		// } else if (correlation.back_strip[0] < 48) {
// 		// 	if (!(
// 		// 		(
// 		// 			correlation.front_energy[0] < 19000 && correlation.back_energy[0] < 28000
// 		// 		)
// 		// 		||
// 		// 		(
// 		// 			correlation.front_energy[0] > 29000 && correlation.back_energy[0] > 42000
// 		// 		)
// 		// 	)) return false;

// 		// } else {

// 		// }

// 	}

// 	return true;
// }


// //-----------------------------------------------------------------------------
// // 									T0D2
// //-----------------------------------------------------------------------------

// T0D2::T0D2(int run)
// : DSSD(run, "t0d2")
// {
// }


// bool T0D2::NormalizeFrontEnergyCheck(
// 	const CorrelatedEvent &correlation,
// 	bool iteration
// )
// {
// 	if
// 	(
// 		correlation.front_energy[0] < 5000
// 		|| correlation.back_energy[0] < 5000
// 	)
// 	{
// 		return false;
// 	}

// 	if (!iteration)
// 	{
// 		if
// 		(
// 			(
// 				correlation.front_energy[0] > 20000
// 				&& correlation.front_energy[0] < 34000
// 			)
// 			||
// 			(
// 				correlation.back_energy[0] > 22000
// 				&& correlation.back_energy[0] < 38000
// 			)
// 		)
// 		{
// 			return false;
// 		}
// 	}
// 	else
// 	{
// 		if
// 		(
// 			abs
// 			(
// 				NormalizeEnergy
// 				(
// 					0,
// 					correlation.front_strip[0],
// 					correlation.front_energy[0]
// 				)
// 				- correlation.back_energy[0]
// 			)
// 			> 1000
// 		)
// 		{
// 			return false;
// 		}
// 	}
// 	return true;
// }


// bool T0D2::NormalizeBackEnergyCheck(
// 	const CorrelatedEvent &correlation,
// 	bool iteration
// )
// {
// 	if
// 	(
// 		correlation.front_energy[0] < 5000
// 		|| correlation.back_energy[0] < 5000
// 	)
// 	{
// 		return false;
// 	}

// 	if (!iteration)
// 	{
// 	}
// 	else
// 	{
// 		if
// 		(
// 			abs(
// 				NormalizeEnergy(
// 					0,
// 					correlation.front_strip[0],
// 					correlation.front_energy[0]
// 				)
// 				- NormalizeEnergy(
// 					1,
// 					correlation.back_strip[0],
// 					correlation.back_energy[0]
// 				)
// 			)
// 			> 1000
// 		)
// 		{
// 			return false;
// 		}
// 	}
// 	return true;
// }


// //-----------------------------------------------------------------------------
// // 									T0D3
// //-----------------------------------------------------------------------------

// T0D3::T0D3(int run)
// : DSSD(run, "t0d3") {
// }


// bool T0D3::NormalizeFrontEnergyCheck
// (
// 	const CorrelatedEvent &correlation,
// 	bool iteration
// )
// {
// 	if
// 	(
// 		correlation.front_energy[0] < 2000
// 		|| correlation.back_energy[0] < 2000
// 	)
// 	{
// 		return false;
// 	}

// 	if (!iteration)
// 	{
// 	}
// 	else
// 	{
// 		if
// 		(
// 			abs
// 			(
// 				NormalizeEnergy
// 				(
// 					0,
// 					correlation.front_strip[0],
// 					correlation.front_energy[0]
// 				)
// 				- correlation.back_energy[0]
// 			)
// 			> 1000
// 		)
// 		{
// 			return false;
// 		}
// 	}
// 	return true;
// }


// bool T0D3::NormalizeBackEnergyCheck
// (
// 	const CorrelatedEvent &correlation,
// 	bool iteration
// )
// {
// 	if
// 	(
// 		correlation.front_energy[0] < 2000
// 		|| correlation.back_energy[0] < 2000
// 	)
// 	{
// 		return false;
// 	}

// 	if (!iteration)
// 	{
// 	}
// 	else
// 	{
// 		if
// 		(
// 			abs(
// 				NormalizeEnergy(
// 					0,
// 					correlation.front_strip[0],
// 					correlation.front_energy[0]
// 				)
// 				- NormalizeEnergy(
// 					1,
// 					correlation.back_strip[0],
// 					correlation.back_energy[0]
// 				)
// 			)
// 			> 1000
// 		)
// 		{
// 			return false;
// 		}
// 	}
// 	return true;
// }




}