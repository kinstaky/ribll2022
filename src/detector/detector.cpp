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
	const std::string &name,
	const std::string &tag
)
: run_(run)
, name_(name)
, tag_(tag) {
}

int Detector::ReadTriggerTimes(std::vector<double> &trigger_times) {
	// clear data
	trigger_times.clear();
	// trigger file name
	TString trigger_file_name;
	trigger_file_name.Form(
		"%s%sxt-map-%s%04d.root",
		kGenerateDataPath,
		kMappingDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
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
	printf("Reading %s trigger events   0%%", tag_.c_str());
	fflush(stdout);
	long long entries = trigger_tree->GetEntries();
	// 1/100 of entry, for showing process
	long long entry100 = entries / 100 + 1;
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


int Detector::MatchTrigger(double, double) {
	// do nothing but report error
	std::cerr << "Error: MatchTrigger is not implemented yet.\n";
	return -1;
}


int Detector::ExtractTrigger(double, double) {
	// do nothing but report error
	std::cerr << "Error: ExtractTrigger is not implemented yet.\n";
	return -1;
}



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










}