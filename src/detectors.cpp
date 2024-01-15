#include "include/detectors.h"

#include <iostream>

#include <TString.h>

#include "include/event/dssd_event.h"
#include "include/event/trigger_event.h"

namespace ribll {

std::shared_ptr<Detector> CreateDetector(
	const std::string &name,
	unsigned int run,
	const std::string &tag
) {
	if (name == "tof" || name == "vtof") {
		return std::make_shared<Tof>(run, name, tag);
	} else if (name == "vt") {
		return std::make_shared<VmeTrigger>(run, tag);
	} else if (name == "xppac") {
		return std::make_shared<Ppac>(run, "xppac", tag);
	} else if (name == "vppac") {
		return std::make_shared<Ppac>(run, "vppac", tag);
	} else if (name == "t0d1") {
		return std::make_shared<T0d1>(run, tag);
	} else if (name == "t0d2") {
		return std::make_shared<T0d2>(run, tag);
	} else if (name == "t0d3") {
		return std::make_shared<T0d3>(run, tag);
	} else if (name.size() == 5 && name.substr(0, 4) == "tafd") {
		unsigned int index = name[4] - '0';
		if (index <= 5) {
			return std::make_shared<Tafd>(run, index, tag);
		}
	} else if (name.size() == 5 && name.substr(0, 4) == "tabd") {
		unsigned int index = name[4] - '0';
		if (index <= 5) {
			return std::make_shared<Tabd>(run, index, tag);
		}
	} else if (name == "tafcsi" || name == "tabcsi") {
		return std::make_shared<CircularCsi>(run, name, tag);
	} else if (name == "t0csi" || name == "t1csi") {
		return std::make_shared<SquareCsi>(run, name, tag);
	} else if (name.size() == 4 && name.substr(0, 3) == "t0s") {
		return std::make_shared<Ssd>(run, name, tag);
	} else if (name == "t1s1") {
		return std::make_shared<Ssd>(run, name, tag);
	}

	std::cerr << "Error: Create detector " << name
		<< " with tag " << tag << " failed.\n";
	return nullptr;
}


std::shared_ptr<Dssd> CreateDssd(
	const std::string &name,
	unsigned int run,
	const std::string &tag
) {
	if (name == "t0d1") {
		return std::make_shared<T0d1>(run, tag);
	} else if (name == "t0d2") {
		return std::make_shared<T0d2>(run, tag);
	} else if (name == "t0d3") {
		return std::make_shared<T0d3>(run, tag);
	} else if (name.size() == 5 && name.substr(0, 4) == "tafd") {
		unsigned int index = name[4] - '0';
		if (index <= 5) {
			return std::make_shared<Tafd>(run, index, tag);
		}
	} else if (name.size() == 4 && name.substr(0, 4) == "tabd") {
		unsigned int index = name[4] - '0';
		if (index <= 5) {
			return std::make_shared<Tabd>(run, index, tag);
		}
	}

	std::cerr << "Error: Create dssd " << name << " failed.\n";
	return nullptr;
}

struct MergeInfo {
	TFile *opf;
	TTree *opt;
	long long entries;
	long long entry;
	std::vector<double> times;
};

int MergeAdssdTrigger(const std::string &trigger_tag, unsigned int run) {
	// VME trigger events
	std::vector<TriggerEvent> vt_events;
	// ADSSD fundamental events
	std::vector<DssdFundamentalEvent> adssd_events[4];

	// merge information of extracted trigger and detector
	std::vector<MergeInfo> merge_info(5);

	// input or output data
	DssdFundamentalEvent fundamental_event;
	TriggerEvent vme_trigger_event;
	double xia_time;

	// read adssd events
	for (int i = 2; i < 6; ++i) {
		// detector input file name
		TString taf_input_file_name;
		taf_input_file_name.Form(
			"%s%stafd%d-fundamental-%stafd%d-%04u.root",
			kGenerateDataPath, kFundamentalDir,
			i, trigger_tag.c_str(), i, run
		);
		// detector input file
		TFile *ipf = new TFile(taf_input_file_name, "read");
		// detector input tree
		TTree *ipt = (TTree*)ipf->Get("tree");
		if (!ipt) {
			std::cerr << "Error: Get tree from "
				<< taf_input_file_name << " failed.\n";
			return -1;
		}

		// XIA trigger file name
		TString xt_file_name;
		xt_file_name.Form(
			"%s%sxt-map-%stafd%d-%04u.root",
			kGenerateDataPath, kMappingDir, trigger_tag.c_str(), i, run
		);
		// add XIA trigger to detector input tree
		ipt->AddFriend("xt=tree", xt_file_name);

		// setup input branches
		ipt->SetBranchAddress("xt.time", &xia_time);
		fundamental_event.SetupInput(ipt);

		MergeInfo &info = merge_info[i-2];
		info.entries = ipt->GetEntries();
		info.entry = 0;
		// read events
		for (long long entry = 0; entry < info.entries; ++entry) {
			ipt->GetEntry(entry);
			info.times.push_back(xia_time);
			adssd_events[i-2].push_back(fundamental_event);
		}

		// close input file
		ipf->Close();
	}

	// read VME trigger events
	// VME trigger input file name
	TString vt_input_file_name;
	vt_input_file_name.Form(
		"%s%svt-fundamental-%svt-%04u.root",
		kGenerateDataPath, kFundamentalDir, trigger_tag.c_str(), run
	);
	TFile *vt_input_file = new TFile(vt_input_file_name, "read");
	TTree *vt_input_tree = (TTree*)vt_input_file->Get("tree");
	if (!vt_input_tree) {
		std::cerr << "Error: Get tree from "
			<< vt_input_file << " failed.\n";
		return -1;
	}
	// XIA trigger input file name
	TString xt_input_file_name;
	xt_input_file_name.Form(
		"%s%sxt-map-%svt-%04u.root",
		kGenerateDataPath, kMappingDir, trigger_tag.c_str(), run
	);
	vt_input_tree->AddFriend("xt=tree", xt_input_file_name);
	// setup input branches
	vt_input_tree->SetBranchAddress("xt.time", &xia_time);
	vme_trigger_event.SetupInput(vt_input_tree);
	// read VME trigger events
	merge_info[4].entries = vt_input_tree->GetEntries();
	merge_info[4].entry = 0;
	for (long long entry = 0; entry < merge_info[4].entries; ++entry) {
		vt_input_tree->GetEntry(entry);
		merge_info[4].times.push_back(xia_time);
		vt_events.push_back(vme_trigger_event);
	}
	// close vt input file
	vt_input_file->Close();


	// setup output file and tree
	for (int i = 2; i < 6; ++i) {
		// detector output file name
		TString taf_output_file_name;
		taf_output_file_name.Form(
			"%s%stafd%d-fundamental-%sta-%04u.root",
			kGenerateDataPath, kFundamentalDir,
			i, trigger_tag.c_str(), run
		);
		// detector output file
		merge_info[i-2].opf = new TFile(taf_output_file_name, "recreate");
		// detector output tree
		merge_info[i-2].opt = new TTree("tree", "fundamental tree");
		// setup output branches
		fundamental_event.SetupOutput(merge_info[i-2].opt);
	}
	// VME output file name
	TString vt_output_file_name;
	vt_output_file_name.Form(
		"%s%svt-fundamental-%sta-%04u.root",
		kGenerateDataPath, kFundamentalDir, trigger_tag.c_str(), run
	);
	merge_info[4].opf = new TFile(vt_output_file_name, "recreate");
	merge_info[4].opt = new TTree("tree", "fundamental tree");
	// setup output branches
	vme_trigger_event.SetupOutput(merge_info[4].opt);

	// setup XIA merged trigger output
	// XIA trigger output file name
	TString xt_output_file_name;
	xt_output_file_name.Form(
		"%s%sxt-map-%sta-%04u.root",
		kGenerateDataPath, kMappingDir, trigger_tag.c_str(), run
	);
	// XIA merged trigger output file
	TFile *xt_output_file = new TFile(xt_output_file_name, "recreate");
	// XIA extracted trigger output tree
	TTree *xt_output_tree = new TTree("tree", "xt tree");
	// setup xt output branch
	xt_output_tree->Branch("time", &xia_time, "t/D");

	// number of finished tree
	size_t finished = 0;
	while (finished < merge_info.size()) {
		double min_time = -1;
		// loop to find the minimum trigger time
		for (const auto &info : merge_info) {
			if (info.entry >= info.entries) continue;
			if (min_time < 0 || min_time > info.times[info.entry]) {
				min_time = info.times[info.entry];
			}
		}
		// loop to fill events
		for (size_t i = 0; i < merge_info.size(); i++) {
			MergeInfo &info = merge_info[i];
			double time = info.times[info.entry];
			if (info.entry >= info.entries || min_time != time) {
				// fill null event
				fundamental_event.Nullify();
				vme_trigger_event.time = -1e5;
				info.opf->cd();
				info.opt->Fill();
			} else {
				// fill valid event
				if (i == 4) {
					vme_trigger_event = vt_events[info.entry];
				} else {
					fundamental_event = adssd_events[i][info.entry];
				}
				info.opt->Fill();
				// move to next entry
				++info.entry;
				if (info.entry == info.entries) {
					++finished;
				}
			}
		}
		// fill output merged XIA trigger
		xia_time = min_time;
		xt_output_tree->Fill();
	}

	// write trees and close files
	xt_output_file->cd();
	xt_output_tree->Write();
	xt_output_file->Close();
	for (auto &info : merge_info) {
		info.opf->cd();
		info.opt->Write();
		info.opf->Close();
	}

	MatchTriggerStatistics statistics(
		run, "ta", trigger_tag, "", 1, 1
	);
	statistics.Write();

	return 0;
}


int MatchWithoutTrigger(const std::string &detector, unsigned int run) {
	const double look_window = 10'000;
	const double window = 100;

	// input map file name
	TString map_file_name;
	map_file_name.Form(
		"%s%s%s-map-%04u.root",
		kGenerateDataPath,
		kMappingDir,
		detector.c_str(),
		run
	);
	// input file
	TFile ipf(map_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< map_file_name << " failed.\n";
		return -1;
	}
	// map event
	DssdMapEvent map_event;
	// setup input branches
	map_event.SetupInput(ipt);

	std::multimap<double, DssdMapEvent> events;

	// read events
	long long entries = ipt->GetEntries();
	long long entry100 = entries / 100;
	// show start
	printf("Reading map events   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get entry
		ipt->GetEntry(entry);
		// record time
		events.insert(std::make_pair(map_event.time, map_event));
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// close file
	ipf.Close();

	// name of output file
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-fundamental-%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		detector.c_str(),
		// tag_.empty() ? "" : (tag_+"-").c_str(),
		run
	);
	// pointer to output file
	TFile opf(output_file_name, "recreate");
	// histogram of look window
	TH1F hist_look_window(
		"ht", "look window", 1000, 0, look_window
	);
	// pointer to output tree
	TTree opt("tree", "fundamental events matched without trigger");
	// ouput event
	DssdFundamentalEvent fundamental_event;
	// setup output branches
	fundamental_event.SetupOutput(&opt);
	// for convenient
	int &fhit = fundamental_event.front_hit;
	int &bhit = fundamental_event.back_hit;

	// show start
	printf("Filling fundamental events   0%%");
	fflush(stdout);
	long long entry = 0;
	for (auto ievent = events.begin(); ievent != events.end(); ++ievent) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		entry++;
		// jump used events
		if (ievent->second.side == 2) continue;
		// initialize
		fhit = 0;
		bhit = 0;
		if (ievent->second.side == 0) {
			fundamental_event.front_strip[fhit] = ievent->second.strip;
			fundamental_event.front_time[fhit] = ievent->second.time;
			fundamental_event.front_energy[fhit] = ievent->second.energy;
			fundamental_event.cfd_flag |=
				ievent->second.cfd_flag ? (1 << fhit) : 0;
			fundamental_event.front_decode_entry[fhit] =
				ievent->second.decode_entry;
			++fhit;
		} else {
			fundamental_event.back_strip[bhit] = ievent->second.strip;
			fundamental_event.back_time[bhit] = ievent->second.time;
			fundamental_event.back_energy[bhit] = ievent->second.energy;
			fundamental_event.cfd_flag |=
				ievent->second.cfd_flag ? (1 << (bhit+8)) : 0;
			fundamental_event.back_decode_entry[bhit] =
				ievent->second.decode_entry;
			++bhit;
		}
		// search correlated events
		for (
			auto jevent = std::next(ievent);
			jevent != events.end();
			++jevent
		) {
			// fill to look window
			hist_look_window.Fill(jevent->first - ievent->first);
			// out of look range
			if (jevent->first >= ievent->first + look_window) break;
			// out of correlated range
			if (jevent->first >= ievent->first + window) continue;
			// jump used events
			if (jevent->second.side == 2) continue;
			if (jevent->second.side == 0) {
				// fill front event
				jevent->second.side = 2;
				if (fhit == 8) continue;
				fundamental_event.front_strip[fhit] = jevent->second.strip;
				fundamental_event.front_time[fhit] = jevent->second.time;
				fundamental_event.front_energy[fhit] = jevent->second.energy;
				fundamental_event.cfd_flag |=
					jevent->second.cfd_flag ? (1 << fhit) : 0;
				fundamental_event.front_decode_entry[fhit] =
					jevent->second.decode_entry;
				++fhit;
			} else {
				// fill back event
				jevent->second.side = 2;
				if (bhit == 8) continue;
				fundamental_event.back_strip[bhit] = jevent->second.strip;
				fundamental_event.back_time[bhit] = jevent->second.time;
				fundamental_event.back_energy[bhit] = jevent->second.energy;
				fundamental_event.cfd_flag |=
					jevent->second.cfd_flag ? (1 << (bhit+8)) : 0;
				fundamental_event.back_decode_entry[bhit] =
					jevent->second.decode_entry;
				++bhit;
			}
		}
		// sort events by energy
		fundamental_event.Sort();
		// fill tree
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save histograms
	hist_look_window.Write();
	// close file
	opf.Close();

	return 0;
}



}	// namespace ribll