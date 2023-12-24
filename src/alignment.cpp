#include "include/alignment.h"

#include <iostream>
#include <iomanip>
#include <algorithm>

#include <TGraph.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>

#include "include/defs.h"
#include "include/statistics/align_statistics.h"

namespace ribll {

Alignment::Alignment(
	unsigned int run,
	size_t group_num,
	double search_window,
	double search_low_bound,
	double search_high_bound
)
: run_(run)
, group_num_(group_num)
, search_window_(search_window)
, search_low_bound_(search_low_bound)
, search_high_bound_(search_high_bound)
, verbose_(false) {
}


Alignment::~Alignment() {
}


int Alignment::ReadXiaTime() {
	// XIA input file name
	TString file_name;
	file_name.Form("%s%svt-map-%04d.root", kGenerateDataPath, kMappingDir, run_);
	// XIA input file
	TFile ipf(file_name, "read");
	// XIA input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: get tree from " << file_name << " failed.\n";
		return -1;
	}
	// set branches
	long long xia_time;
	ipt->SetBranchAddress("timestamp", &xia_time);

	// initialize array
	xia_times_.clear();

	// show start
	printf("reading xia events   0%%");
	fflush(stdout);
	// 1/100 of total entries
	long long entry100 = ipt->GetEntries() / 100 + 1;
	for (long long entry = 0; entry < ipt->GetEntries(); ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);
		xia_times_.push_back(xia_time);
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// sort the timestamp from xia
	std::sort(xia_times_.begin(), xia_times_.end());
	// close input file
	ipf.Close();
	return 0;
}


int Alignment::ReadVmeTime() {
	// VME input file name
	TString file_name;
	file_name.Form("%s%s%04d.root", kCrate3Path, kCrate3FileName, run_);
	// VME event input file
	TFile ipf(file_name, "read");
	// VME input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: get tree from " << file_name << " failed.\n";
		return -1;
	}
	// set branches
	unsigned long long sdc[32];
	ipt->SetBranchAddress("sdc", sdc);

	// initialize
	vme_times_.clear();
	// offset cause by bit flip
	unsigned long long bit_flip_offset = 0;
	// last timestamp for checking bit flip
	unsigned long long last_timestamp = 0;
	// 1/100 of total entries
	long long entry100 = ipt->GetEntries() / 100 + 1;
	// show start
	printf("reading vme events   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < ipt->GetEntries(); ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);
		// calculate timestamp
		unsigned long long timestamp =
			(sdc[kScalerIndex] + bit_flip_offset) * kScalerPeriod;
		// check sdc over range
		if (timestamp < last_timestamp) {
			bit_flip_offset += 1llu << 32;
			timestamp += (1llu << 32) * kScalerPeriod;
		}
		last_timestamp = timestamp;

		// record it
		vme_times_.push_back(timestamp);
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// close input file
	ipf.Close();
	return 0;
}

struct Edge {
	// left edge
	double left;
	// right edge
	double right;
};


Edge EdgeDetect(TH1F *hist) {
	const double threshold = 3;
	const int check_num = 4;
	// detect rise edge
	bool detect_rise = false;
	// detect fall edge
	bool detect_fall = false;
	// return result of left edge and right edge
	Edge result{0.0, 0.0};
	// bins of histogram
	int bins = hist->GetNbinsX();
	for (int i = 1; i <= bins-check_num+1; ++i) {
		int over_threshold = 0;
		for (int j = 0; j < check_num; ++j) {
			if (hist->GetBinContent(i+j) > threshold) {
				++over_threshold;
			}
		}
		if (!detect_rise) {
			if (over_threshold == check_num) {
				detect_rise = true;
				result.left = hist->GetBinCenter(i);
				result.right = result.left;
			}
		} else if (over_threshold == 0) {
			detect_fall = true;
			result.right = hist->GetBinCenter(i);
			break;
		}
	}
	if (!detect_fall) result.right = hist->GetBinCenter(bins-check_num+1);
	return result;
}

void Alignment::GroupAlignment() {
	// offset-group graph in first, second. third match
	TGraph offset_vs_group[3];
	// window-group graph, in first, second, third match
	TGraph window_vs_group[3];
	// output tree
	TTree *opt = new TTree("tree", "alignment result");
	// output match VME trigger time recorded in XIA
	long long xia_time;
	// output vme trigger time recorded in VME
	long long vme_time;
	// setup branches
	opt->Branch("xia_time", &xia_time, "xt/L");
	opt->Branch("vme_time", &vme_time, "vt/L");

	// initialize statistics
	statistics_ = AlignStatistics(
		run_, xia_times_.size(), vme_times_.size(),  group_num_
	);

	// group size
	size_t group_size = vme_times_.size() / group_num_ + 1;
	// lower bound of offset in each group, set general for the first group
	double lower_bound = search_low_bound_;
	// upper bound of offset in each group, set general for the first group
	double upper_bound = search_high_bound_;
	for (int group = 0; group < group_num_; ++group) {
		if (verbose_) {
			std::cout << "----------------------------------------\n"
				<< "group " << group << "\n";
		}

		// start VME event index of this group
		size_t group_start_index = group_size * group;
		// end VME event index of this group
		size_t group_end_index = group_size * (group + 1);
		group_end_index = group_end_index < vme_times_.size()
			? group_end_index : vme_times_.size();

		// maximum match rate in all offsets
		double max_match_rate = 0;
		// offset of max match rate
		double max_match_offset = 0;
		// graph of match rate VS offset
		TGraph match_rate_vs_offset;
		// loop to search the for the offset with fixed steps
		for (
			double offset = lower_bound;
			offset <= upper_bound;
			offset += search_window_
		) {
			// number of matched events
			long long match_events = 0;
			// first time match to search for the max match rate
			for (
				size_t vme_entry = group_start_index;
				vme_entry < group_end_index;
				++vme_entry
			) {
				// break the loop if reach the end of array
				if (vme_entry >= vme_times_.size()) break;
				// VME time
				long long vme_time = vme_times_[vme_entry];
				// loop XIA events to find match events
				for (
					auto xia_iter = std::lower_bound(
						xia_times_.begin(),
						xia_times_.end(),
						vme_time + offset
					);
					xia_iter != xia_times_.end();
					++xia_iter
				) {
					if (*xia_iter > vme_time + offset + search_window_) break;
					++match_events;
				}
			}
			// record match rate in graph
			double match_rate = double(match_events) / double(group_size);
			match_rate_vs_offset.AddPoint(offset, match_rate);
			if (max_match_rate < match_rate) {
				max_match_rate = match_rate;
				max_match_offset = offset;
			}
			// show match rate of each seconds
			if (verbose_ && (long long)offset % 1'000'000'000 == 0) {
				std::cout << "  " << std::setprecision(10)
					<< offset / 1'000'000'000.0 << " s"
					<< "  max match rate " << max_match_rate
					<< "  offset " << max_match_offset << "\n";
			}
		}
		// save rate-offset graph
		match_rate_vs_offset.Write(TString::Format("rate%d", group));

		// current offset
		double offset = max_match_offset;
		// current window
		double window = search_window_;
		// total number of match events
		long long match_events = 0;
		// oversize events
		long long oversize = 0;
		// write to graph
		offset_vs_group[0].AddPoint(group, offset);
		window_vs_group[0].AddPoint(group, window);

		// time window for first match
		TH1F first_match_window(
			TString::Format("ht%d_1", group), "window", 1000, 0, window
		);
		// first match out of the loop to calculate the average and variance
		for (
			size_t vme_entry = group_start_index;
			vme_entry < group_end_index;
			++vme_entry
		) {
			// break the loop if reach the end of array
			if (vme_entry >= vme_times_.size()) break;
			// VME time
			long long vme_time = vme_times_[vme_entry];
			// match count of current VME event
			int match_count = 0;
			// loop XIA events to find match events
			for (
				auto xia_iter = std::lower_bound(
					xia_times_.begin(),
					xia_times_.end(),
					vme_time + offset
				);
				xia_iter != xia_times_.end();
				++xia_iter
			) {
				if (*xia_iter > vme_time + offset + window) break;
				double difference = double(*xia_iter - vme_time - offset);
				first_match_window.Fill(difference);
				++match_events;
				++match_count;
			}
			if (match_count == 1) {
				++statistics_.first_align_events;
			} else {
				++oversize;
			}
		}
		// detect edge
		Edge edge = EdgeDetect(&first_match_window);
		// save time window
		first_match_window.Write();
		// show the first match result
		if (verbose_) {
			std::cout << "group " << group
				<< ", first match offset " << offset
				<< ", window " << window
				<< ", rate " << double(match_events) / group_size
				<< ", oversize " << oversize << "\n"
				<< "edge left " << edge.left
				<< ", edge right " << edge.right << "\n";
		}

		// the second match to adjust the offset and window
		// initialize
		offset += (edge.left + edge.right) / 2.0;
		window = (edge.right - edge.left) * 2.0;
		match_events = 0;
		oversize = 0;
		// write to graph
		offset_vs_group[1].AddPoint(group, offset);
		window_vs_group[1].AddPoint(group, window);
		// time window in the second match
		TH1F second_match_window(
			TString::Format("ht%d_2", group), "window", 50, -window, window
		);
		// the second match
		for (
			size_t vme_entry = group_start_index;
			vme_entry < group_end_index;
			++vme_entry
		) {
			// break the loop if reach the end of array
			if (vme_entry >= vme_times_.size()) break;
			// VME time
			long long vme_time = vme_times_[vme_entry];
			// match count of current VME event
			int match_count = 0;
			// loop XIA events to find match events
			for (
				auto xia_iter = std::lower_bound(
					xia_times_.begin(),
					xia_times_.end(),
					vme_time + offset - window
				);
				xia_iter != xia_times_.end();
				++xia_iter
			) {
				if (*xia_iter > vme_time + offset + window) break;
				double difference = double(*xia_iter - vme_time - offset);
				second_match_window.Fill(difference);
				++match_events;
				++match_count;
			}
			if (match_count == 1) {
				++statistics_.second_align_events;
			} else {
				++oversize;
			}
		}
		// detect edge
		edge = EdgeDetect(&second_match_window);
		// save time window
		second_match_window.Write();

		// show the second match result
		if (verbose_) {
			std::cout
				<< "group " << group
				<< ", second match offset " << offset
				<< ", window " << window
				<< ", rate " << double(match_events) / group_size
				<< ", oversize " << oversize << "\n"
				<< "edge left " << edge.left
				<< ", edge right " << edge.right << "\n";
		}

		// the third match to fill the tree
		// initialize
		offset += (edge.left + edge.right) / 2.0;
		window = (edge.right - edge.left) / 1.8;
		match_events = 0;
		oversize = 0;
		// write to graph
		offset_vs_group[2].AddPoint(group, offset);
		window_vs_group[2].AddPoint(group, window);
		// time window in the third match
		TH1F third_match_window(
			TString::Format("ht%d_3", group), "window", 30, -window, window
		);
		// the third match to fill the tree
		for (
			size_t vme_entry = group_start_index;
			vme_entry < group_end_index;
			++vme_entry
		) {
			// break the loop if reach the end of array
			if (vme_entry >= vme_times_.size()) break;
			// VME time
			vme_time = vme_times_[vme_entry];
			// match count of current VME event
			int match_count = 0;
			// initialize XIA time
			xia_time = -1.0;
			// loop XIA events to find match events
			for (
				auto xia_iter = std::lower_bound(
					xia_times_.begin(),
					xia_times_.end(),
					vme_time + offset - window
				);
				xia_iter != xia_times_.end();
				++xia_iter
			) {
				if (*xia_iter > vme_time + offset + window) break;
				double difference = double(*xia_iter - vme_time - offset);
				third_match_window.Fill(difference);
				++match_events;
				++match_count;
				xia_time = *xia_iter;
			}
			if (match_count == 1) {
				++statistics_.third_align_events;
				++statistics_.align_events;
			} else if (match_count > 1) {
				++oversize;
				++statistics_.oversize_events;
			}
			opt->Fill();
		}
		// save time window
		third_match_window.Write();
		// show the third match result
		if (verbose_) {
			std::cout
				<< "group " << group
				<< ", third match offset " << offset
				<< ", window " << window
				<< ", rate " << double(match_events) / group_size
				<< ", oversize " << oversize << "\n"
				<< "----------------------------------------\n";
		}

		// upate lower bound and upper bound
		lower_bound = max_match_offset - 2 * search_window_;
		upper_bound = max_match_offset + 2 * search_window_;
	}
	// save tree
	opt->Write();
	// save offsets and windows
	for (int i = 0; i < 3; ++i) {
		offset_vs_group[i].Write(TString::Format("goffset%d", i));
		window_vs_group[i].Write(TString::Format("gwindow%d", i));
	}
	// save and print statistics
	statistics_.Write();
	statistics_.Print();
	return;
}


int GetGdcOffset(unsigned int run) {
	// input VME file name
	TString vme_file_name;
	vme_file_name.Form(
		"%s%s%04d.root",
		kCrate3Path, kCrate3FileName, run
	);
	// input VME file
	TFile vme_file(vme_file_name, "read");
	// input VME tree
	TTree *vme_tree = (TTree*)vme_file.Get("tree");
	if (!vme_tree) {
		std::cerr << "Error: Get tree from "
			<< vme_file_name << " failed.\n";
		return -1;
	}
	// input data
	int gmulti[2][128];
	int madc[2][32];
	// setup input branches
	vme_tree->SetBranchAddress("gmulti", gmulti);
	vme_tree->SetBranchAddress("madc", madc);
	std::vector<bool> time;
	std::vector<bool> energy;

	// total number of entries
	long long entries = vme_tree->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Reading values   0%%");
	fflush(stdout);
	// read events
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		vme_tree->GetEntry(entry);

		// read values into arrays
		time.push_back(gmulti[1][0] > 0 ? true : false);
		energy.push_back(madc[0][16] > 200 ? true  : false);
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// close input file
	vme_file.Close();

	std::vector<int> offset_match;
	int size = energy.size();
	for (int offset = 0; offset < 1500; ++offset) {
		int match = 0;
		for (int i = 0; i < size; ++i) {
			if (i + offset < 0) continue;
			if (i + offset >= size) continue;
			if (energy[i] && time[i+offset]) ++match;
		}
		offset_match.push_back(match);
	}

	// find max match in GDC
	// max offset
	int max_offset = 0;
	// max match number
	int max_match = offset_match[0];
	for (size_t i = 1; i < offset_match.size(); ++i) {
		if (offset_match[i] > max_match) {
			max_match = offset_match[i];
			max_offset = i;
		}
	}
	std::cout << "GDC 1: " << max_offset << "\n";
	return max_offset;
}


int Alignment::AlignGdc() {
	// this run input file name
	TString first_input_file_name;
	first_input_file_name.Form(
		"%s%s%04u.root",
		kCrate3Path, kCrate3FileName, run_
	);
	// this run input file
	TFile *first_input_file = new TFile(first_input_file_name, "read");
	// this run input tree
	TTree *first_tree = (TTree*)first_input_file->Get("tree");
	if (!first_tree) {
		std::cerr << "Error: Get tree from "
			<< first_input_file_name << " failed.\n";
		return -1;
	}
	// input gmulti
	int gmulti[2][128];
	// input gdc
	int gdc[2][128][5];
	// setup input branches
	first_tree->SetBranchAddress("gmulti", gmulti);
	first_tree->SetBranchAddress("gdc", gdc);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%salign-gdc-%04u.root",
		kGenerateDataPath, kAlignDir, run_
	);
	// output file
	TFile *opf = new TFile(output_file_name, "recreate");
	// output tree
	TTree *opt = new TTree("tree", "align gdc");
	// output gmulti
	int align_gmulti[128];
	// output gdc
	int align_gdc[128][5];
	// setup output branches
	opt->Branch("gmulti", align_gmulti, "gmulti[128]/I");
	opt->Branch("gdc", align_gdc, "gdc[128][5]/I");

	// total entries
	long long entries = first_tree->GetEntries();

	long long offset = GetGdcOffset(run_);
	if (offset < 0) {
		std::cerr << "Error: Searching gdc 1 offset failed.\n";
		return -1;
	}

	// first loop to read from first input file (this run)
	for (long long entry = offset; entry < entries; ++entry) {
		first_tree->GetEntry(entry);
		for (int i = 0; i < 128; ++i) {
			align_gmulti[i] = gmulti[1][i];
			for (int j = 0; j < 5; ++j) {
				align_gdc[i][j] = gdc[1][i][j];
			}
		}
		opt->Fill();
	}

	// close first file
	first_input_file->Close();


	// next run input file name
	TString second_input_file_name;
	second_input_file_name.Form(
		"%s%s%04u.root",
		kCrate3Path, kCrate3FileName, run_+1
	);
	// this run input file
	TFile *second_input_file = new TFile(second_input_file_name, "read");
	// this run input tree
	TTree *second_tree = (TTree*)second_input_file->Get("tree");
	if (!second_tree) {
		std::cerr << "Error: Get tree from "
			<< second_input_file_name << " failed.\n";
		return -1;
	}
	// setup input branches
	second_tree->SetBranchAddress("gmulti", gmulti);
	second_tree->SetBranchAddress("gdc", gdc);

	// second loop to read from second input file (next run)
	for (long long entry = 0; entry < offset; ++entry) {
		second_tree->GetEntry(entry);
		for (int i = 0; i < 128; ++i) {
			align_gmulti[i] = gmulti[1][i];
			for (int j = 0; j < 5; ++j) {
				align_gdc[i][j] = gdc[1][i][j];
			}
		}
		opt->Fill();
	}

	// close first file
	second_input_file->Close();

	// save align time and close files
	opf->cd();
	opt->Write();
	opf->Close();

	return 0;
}


int Alignment::Align() {
	// read xia and vme timestamp
	if (ReadXiaTime()) {
		std::cerr << "Error: read xia time failed.\n";
		return -1;
	}
	if (ReadVmeTime()) {
		std::cerr << "Error: read vme time failed.\n";
		return -1;
	}

	// setup output alignment cache file
	TString output_file_name;
	output_file_name.Form(
		"%s%salign-%04d.root",
		kGenerateDataPath, kAlignDir, run_
	);
	// output file
	TFile output_file(output_file_name, "recreate");
	GroupAlignment();
	output_file.Close();

	// align gdc
	if (AlignGdc()) {
		return -1;
	}
	return 0;
}


}