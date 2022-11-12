#include "include/detector/detector.h"

#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

#include "include/defs.h"

const unsigned int look_window = 10'000;

namespace ribll {

Detector::Detector(
	int run,
	const std::string &name,
	unsigned int single_side_window,
	unsigned int double_sides_window
)
: run_(run)
, name_(name)
, single_side_window_(single_side_window)
, double_sides_window_(double_sides_window) {

	single_side_file_ = nullptr;
	correlated_file_ = nullptr;
	merged_file_ = nullptr;
}


Detector::~Detector() {
	if (single_side_file_) {
		single_side_file_->Close();
	}
	if (correlated_file_) {
		correlated_file_->Close();
	}
	if (merged_file_) {
		merged_file_->Close();
	}
}


int Detector::SetupSingleSideOutput() {
	// single side event file name
	TString single_side_file_name;
	single_side_file_name.Form(
		"%s%s%s-single-side-%04d.root",
		kGenerateDataPath, kSingleSideDir, name_.c_str(), run_
	);
	// setup single side file
	single_side_file_ = new TFile(single_side_file_name, "recreate");
	// setup single side tree
	single_side_tree_ = new TTree(
		"tree", TString::Format("residual single side events of %s", name_.c_str())
	);
	// setup single side tree branches
	single_side_tree_->Branch("index", &single_side_.index, "index/s");
	single_side_tree_->Branch("side", &single_side_.side, "side/s");
	single_side_tree_->Branch("hit", &single_side_.hit, "hit/s");
	single_side_tree_->Branch("strip", single_side_.strip, "s[hit]/s");
	single_side_tree_->Branch("time", single_side_.time, "t[hit]/D");
	single_side_tree_->Branch("energy", single_side_.energy, "e[hit]/D");
	
	return 0;
}


int Detector::SetupCorrelatedOutput() {
	// setup output correlated file
	TString correlated_file_name;
	correlated_file_name.Form(
		"%s%s%s-corr-%04d.root",
		kGenerateDataPath, kCorrelationDir, name_.c_str(), run_
	);
	correlated_file_ = new TFile(correlated_file_name, "recreate");
	// setup output tree
	correlated_tree_ = new TTree(
		"tree", TString::Format("correlated tree of %s", name_.c_str())
	);
	// setup output data
	SetupCorrelatedTree();
	return 0;
}


int Detector::ReadEvents() {
	// setup input file
	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-map-%04d.root", kGenerateDataPath, kMappingDir, name_.c_str(), run_
	);
	TFile *input_file = new TFile(input_file_name, "read");
	TTree *ipt = (TTree*)input_file->Get("tree");
	if (!ipt) {
		std::cerr << "Error: get tree from " << input_file_name << " failed.\n";
		return -1;
	}
	// input data and branches
	long long timestamp;
	unsigned short index, side, strip;
	double time, energy;
	// branches
	ipt->SetBranchAddress("timestamp", &timestamp);
	ipt->SetBranchAddress("index", &index);
	ipt->SetBranchAddress("side", &side);
	ipt->SetBranchAddress("strip", &strip);
	ipt->SetBranchAddress("time", &time);
	ipt->SetBranchAddress("energy", &energy);

	if (SetupSingleSideOutput()) {
		std::cerr << "Error: Setup single side output failed.\n";
		return -1;
	}
	TH1F *single_side_look_window = new TH1F("ht", "look window", 1000, -look_window, look_window);
	TH1F *x_side_look_window = new TH1F("hxt", "look window", 1000, -look_window, look_window);
	TH1F *y_side_look_window = new TH1F("hyt", "look window", 1000, -look_window, look_window);

	// initialize
	events_.clear();

	// show process
	printf("Reading events   0%%");
	fflush(stdout);
	// 1/100 of entry
	long long entry100 = ipt->GetEntries() / 100;

	// to show hit rate at the end of this function
	// total event number of x side
	long long x_total_event = 0;
	// total event number of y side
	long long y_total_event = 0;
	// to show x side single rate
	long long x_single_event = 0;
	// to show y side single rate
	long long y_single_event = 0;
	// to show x sdie double rate
	long long x_double_event = 0;
	// to show y sdie double rate
	long long y_double_event = 0;
	// to show x side triple rate
	long long x_triple_event = 0;
	// to show y side triple rate
	long long y_triple_event = 0;

	for (long long entry = 0; entry < 1'000'000; ++entry) {
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);
	
		// iterator points to the element to erase, to make the key be the 
		// timestmap of the event with largest energy
		auto erase_iter = events_.end();
		// search if in the time-window of exist event
		bool found = false;
		for (
			auto iter = events_.lower_bound(timestamp - look_window);
			iter != events_.lower_bound(timestamp + look_window);
			++iter
		) {
			// ignore other detector event (different index)
			if (iter->second.index != index) continue;
			// ignore other side event
			if (iter->second.side != side) continue;
			// fill to look window
			single_side_look_window->Fill(iter->first - timestamp);
			if (!side) {
				x_side_look_window->Fill(iter->first - timestamp);
			} else {
				y_side_look_window->Fill(iter->first - timestamp);
			}
			// continue if out of single side window
			if (iter->first <= timestamp - single_side_window_) continue;
			if (iter->first >= timestamp + single_side_window_) continue;
			// jump if has found correlated event
			if (found) continue;
			// found correlated event
			found = true;
			if (iter->second.hit < 8) {
				iter->second.strip[iter->second.hit] = strip;
				iter->second.time[iter->second.hit] = time;
				iter->second.energy[iter->second.hit] = energy;
				// use the largest energy event's timestamp
				bool energy_max = true;
				for (size_t i = 0; i < iter->second.hit; ++i) {
					if (energy < iter->second.energy[i]) {
						energy_max = false;
						break;
					}
				}
				// new event has the largest erergy
				// use its timestamp to represent the whole event
				if (energy_max) {
					erase_iter = iter;
				}
			}
			++(iter->second.hit);
			// record single, double and triple event counts
			switch (iter->second.hit) {
				case 2:
					if (iter->second.side == 0) {
						// x side
						--x_single_event;
						++x_double_event;
					} else {
						// y side
						--y_single_event;
						++y_double_event;
					}
					break;
				case 3:
					if (iter->second.side == 0) {
						// x side
						--x_double_event;
						++x_triple_event;
					} else {
						// y side
						--y_double_event;
						++y_triple_event;
					}
					break;
				default:
					// over 3
					if (iter->second.side == 0) {
						// x side
						--x_triple_event;
					} else {
						// y side
						--y_triple_event;
					}
			}
		}
		if (found && erase_iter != events_.end()) {
			SingleSideEvent event = erase_iter->second;
			events_.erase(erase_iter);
			events_.insert(std::make_pair(timestamp, event));
		}
		if (!found) {
			// not found, insert new element to map
			SingleSideEvent event {
				false,
				index,
				side,
				1,
				{strip, 0, 0, 0, 0, 0, 0, 0},
				{time, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
				{energy, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
			};
			events_.insert(std::make_pair(timestamp, event));
			// increase single event
			if (side == 0) {
				// x side
				++x_single_event;
				++x_total_event;
			} else {
				// y side
				++y_single_event;
				++y_total_event;
			}
		}
	}
	printf("\b\b\b\b100%%\n");

	// write look window
	single_side_file_->cd();
	single_side_look_window->Write();
	x_side_look_window->Write();
	y_side_look_window->Write();
	// close input file
	input_file->Close();

	// show rate
	std::cout << "x side total " << x_total_event << "\n"
		<< "  single rate " << x_single_event << " " << double(x_single_event) / x_total_event << "\n"
		<< "  double rate " << x_double_event << " " << double(x_double_event) / x_total_event << "\n"
		<< "  triple rate " << x_triple_event << " " << double(x_triple_event) / x_total_event << "\n"
		<< "y side total " << y_total_event << "\n"
		<< "  single rate " << y_single_event << " " << double(y_single_event) / y_total_event << "\n"
		<< "  double rate " << y_double_event << " " << double(y_double_event) / y_total_event << "\n"
		<< "  triple rate " << y_triple_event << " " << double(y_triple_event) / y_total_event << "\n";

	return 0;
}


int Detector::Correlate() {
	// read detector events
	if (ReadEvents()) {
		std::cerr << "Error: Read " << name_ << " events failed.\n";
		return -1;
	}

	if (SetupCorrelatedOutput()) {
		std::cerr << "Error: Setup correlated output file or tree failed.\n";
		return -1;
	}
	correlated_file_->cd();
	TH1F *double_side_look_window = new TH1F("ht", "look window", 1000, 0, look_window);

	
	// show process
	printf("Correlating %s   0%%", name_.c_str());
	fflush(stdout);
	// 1/100 of entry
	long long entry100 = events_.size() / 100;
	// correlating entry 
	long long entry = 0;
	
	// for calculating the correlated rate
	// correlated rate of events
	long long correlated_events = 0;
	// over correlated events, one to serveral condition
	long long over_correlated_events = 0;
	
	for (auto ievent = events_.begin(); ievent != events_.end(); ++ievent) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		++entry;

		// jump correlated events
		if (ievent->second.used) continue;
		// jump if hit over 8
		if (ievent->second.hit >= 8) continue;


		// found correlated event on the other side
		bool found = false;
		for (
			auto jevent = events_.lower_bound(ievent->first);
			jevent != events_.upper_bound(ievent->first + look_window);
			++jevent
		) {
			// ignore itself
			if (jevent == ievent) continue;
			// fill to look window
			double_side_look_window->Fill(jevent->first - ievent->first);
			
			// jump if out of search window
			if (jevent->first - ievent->first > double_sides_window_) continue;
			// jump correlated event
			if (jevent->second.used) continue;
			// jump other detector
			if (jevent->second.index != ievent->second.index) continue;
			// jump same side event
			if (jevent->second.side == ievent->second.side) continue;
			// jump if hit over 8
			if (jevent->second.hit >= 8) continue;

			if (!found) {
				found = true;
				jevent->second.used = true;
				// fill to correlation_
				// x side event
				SingleSideEvent *x_event = nullptr;
				// y side event
				SingleSideEvent *y_event = nullptr;
				if (ievent->second.side == 0) {
					x_event = &(ievent->second);
					y_event = &(jevent->second);
					correlation_.timestamp = ievent->first;
				} else {
					x_event = &(jevent->second);
					y_event = &(ievent->second);
					correlation_.timestamp = jevent->first;
				}

				correlation_.index = x_event->index;
				// sort event in single side by energy
				correlation_.x_hit = x_event->hit;
				for (unsigned short i = 0; i < x_event->hit; ++i) {
					correlation_.x_strip[i] = x_event->strip[i];
					correlation_.x_time[i] = x_event->time[i];
					correlation_.x_energy[i] = x_event->energy[i];
				}
				correlation_.y_hit = y_event->hit;
				for (unsigned short i = 0; i < y_event->hit; ++i) {
					correlation_.y_strip[i] = y_event->strip[i];
					correlation_.y_time[i] = y_event->time[i];
					correlation_.y_energy[i] = y_event->energy[i];
				}
			} else {
				++over_correlated_events;
			}

		}

		if (found) {
			ievent->second.used = true;
			correlated_tree_->Fill();
			correlated_events += 2;
		}

	}
	// show finished
	printf("\b\b\b\b100%%\n");

	// write and close output file
	correlated_file_->cd();
	double_side_look_window->Write();
	correlated_tree_->Write();
	correlated_file_->Close();

	// show correlated rate
	std::cout << "correlated rate " << correlated_events << " / " << events_.size()
		<< " " << double(correlated_events) / events_.size() << "\n"
		<< "over correlated " << over_correlated_events << "\n";

	// fill residual single side events
	// show process
	printf("Filling residual single side events   0%%");
	fflush(stdout);
	// filtering index
	size_t filtering_index = 0;
	// 1/100 of events size
	size_t filtering_index_100 = events_.size() / 100;
	single_side_file_->cd();
	for (const auto &[key, value] : events_) {
		// show process
		if (filtering_index % filtering_index_100 == 0) {
			printf("\b\b\b\b%3lu%%", filtering_index / filtering_index_100);
			fflush(stdout);
		}
		++filtering_index;

		if (value.used) continue;
		single_side_ = value;

		if (single_side_.hit > 8) continue;
		// bubble sort by strip
		bool exchange = true;
		while (exchange) {
			exchange = false;
			for (unsigned short i = 0; i < single_side_.hit-1; ++i) {
				if (single_side_.strip[i] > single_side_.strip[i+1]) {
					// exchange
					exchange = true;
					// exchange strip
					unsigned short tmp_strip = single_side_.strip[i];
					single_side_.strip[i] = single_side_.strip[i+1];
					single_side_.strip[i+1] = tmp_strip;
					// excahange time
					double tmp = single_side_.time[i];
					single_side_.time[i] = single_side_.time[i+1];
					single_side_.time[i+1] = tmp;
					// exchange energy
					tmp = single_side_.energy[i];
					single_side_.energy[i] = single_side_.energy[i+1];
					single_side_.energy[i+1] = tmp;
				}
			}
		}
		
		single_side_tree_->Fill();
	}
	// show process
	printf("\b\b\b\b100%%\n");
	// write tree
	single_side_tree_->Write();
	return 0;
}

//-----------------------------------------------------------------------------
// 										DSSD
//-----------------------------------------------------------------------------

DSSD::DSSD(
	int run,
	const std::string &name,
	unsigned int single_side_window,
	unsigned int double_sides_window
)
:Detector(run, name, single_side_window, double_sides_window) {
}


void DSSD::SetupCorrelatedTree() {
	correlated_tree_->Branch("timestamp", &correlation_.timestamp, "ts/L");
	correlated_tree_->Branch("xhit", &correlation_.x_hit, "xhit/s");
	correlated_tree_->Branch("yhit", &correlation_.y_hit, "yhit/s");
	correlated_tree_->Branch("xstrip", correlation_.x_strip, "xs[xhit]/s");
	correlated_tree_->Branch("ystrip", correlation_.y_strip, "ys[yhit]/s");
	correlated_tree_->Branch("xtime", correlation_.x_time, "xt[xhit]/D");
	correlated_tree_->Branch("ytime", correlation_.y_time, "yt[yhit]/D");
	correlated_tree_->Branch("xenergy", correlation_.x_energy, "xe[xhit]/D");
	correlated_tree_->Branch("yenergy", correlation_.y_energy, "ye[yhit]/D");
}


}