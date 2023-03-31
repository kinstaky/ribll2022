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
	// set xia input file name
	TString file_name;
	file_name.Form("%s%svt-map-%04d.root", kGenerateDataPath, kMappingDir, run_);
	// read root file
	TFile *ipf = new TFile(file_name, "read");
	if (!ipf) {
		std::cerr << "Error: open file " << file_name << " failed.\n";
		return -1;
	}
	// get tree
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: get tree from " << file_name << " failed.\n";
		ipf->Close();
		return -1;
	}
	// set branches
	long long xia_time;
	ipt->SetBranchAddress("timestamp", &xia_time);

	// show process
	printf("reading xia events   0%%");
	fflush(stdout);
	Long64_t nentry100 = ipt->GetEntries() / 100 + 1;
	xia_times_.clear();
	for (Long64_t entry = 0; entry < ipt->GetEntries(); ++entry) {
		if (entry % nentry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / nentry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);
		xia_times_.push_back(xia_time);

	}
	printf("\b\b\b\b100%%\n");

	// sort the timestamp from xia
	std::sort(xia_times_.begin(), xia_times_.end());

	ipf->Close();
	return 0;
}


int Alignment::ReadVmeTime() {
	// set vme input file name
	TString file_name;
	file_name.Form("%s%s%04d.root", kCrate3Path, kCrate3FileName, run_);
	// raed root file
	TFile *ipf = new TFile(file_name, "read");
	if (!ipf) {
		std::cerr << "Error: open file " << file_name << " failed.\n";
		return -1;
	}
	// get tree
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: get tree from " << file_name << " failed.\n";
		ipf->Close();
		return -1;
	}
	// set branches
	ULong64_t sdc[32];
	ipt->SetBranchAddress("sdc", sdc);

	// show process
	printf("reading vme events   0%%");
	fflush(stdout);
	Long64_t nentry100 = ipt->GetEntries() / 100 + 1;
	vme_times_.clear();
	ULong64_t bit_flip_offset = 0;
	ULong64_t last_timestamp = 0;
	for (Long64_t entry = 0; entry < ipt->GetEntries(); ++entry) {
		if (entry % nentry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / nentry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);
		// calculate timestamp
		ULong64_t timestamp = (sdc[kScalerIndex] + bit_flip_offset) * kScalerPeriod;
		// check sdc over range
		if (timestamp < last_timestamp) {
			bit_flip_offset += 1llu << 32;
			timestamp += (1llu << 32) * kScalerPeriod;
		}
		last_timestamp = timestamp;

		// record it
		vme_times_.push_back(timestamp);

	}
	printf("\b\b\b\b100%%\n");

	ipf->Close();
	return 0;
}


TGraph* Alignment::GroupAlignment() {
	// graph for xia-vme calibration
	TGraph *result = new TGraph;
	// array of group hit
	TGraph **offset_hit = new TGraph*[group_num_];
	// array of group average hit
	TGraph **offset_average_hit = new TGraph*[group_num_];
	// hit of different group
	TGraph *groups_hit = new TGraph;
	// time window of different group
	TH1F **group_window = new TH1F*[group_num_];
	// offsets of different group
	std::vector<double> group_offset;

	// first group loop to find the group offset
	size_t group_size = vme_times_.size() / group_num_ + 1;
	double last_offset = 0;
	for (int g = 0; g < group_num_; ++g) {
		if (verbose_) {
			std::cout << "group " << g << "\n";
		}
		// first children loop that finding the offset
		offset_hit[g] = new TGraph;
		offset_average_hit[g] = new TGraph;
		double max_hit = 0;
		double max_hit_offset = 0;
		double low_bound = g < 3 ? search_low_bound_ : last_offset - 5 * search_window_;
		double high_bound = g < 3 ? search_high_bound_ : last_offset + 5 * search_window_;
		int group_point = 0;
		for (double offset = low_bound; offset <= high_bound; offset += search_window_) {
			int ievent_count = 0;
			Double_t hit = 0.0;

			for (
				size_t vme_entry = group_size * g;
				vme_entry < group_size * (g + 1) && vme_entry < vme_times_.size();
				++vme_entry
			) {
				Long64_t vme_time = vme_times_[vme_entry];
				++ievent_count;
				for (
					auto xia_iter = std::lower_bound(
						xia_times_.begin(),
						xia_times_.end(),
						vme_time - search_window_ + offset
					);
					xia_iter != xia_times_.end();
					++xia_iter
				) {
					if (*xia_iter > vme_time + search_window_ + offset) break;
					++hit;
				}
			}

			offset_hit[g]->SetPoint(group_point, offset, hit);
			Double_t average_hit = hit / ievent_count;
			if (max_hit < average_hit) {
				max_hit = average_hit;
				max_hit_offset = offset;
			}
			offset_average_hit[g]->SetPoint(group_point, offset, average_hit);
			++group_point;


			if (verbose_ && (Long64_t)offset % 1000000000 == 0) {
				std::cout << "  " << std::setprecision(10)
					<< offset / 1000000000.0 << " s  "
					<< hit << "  " << average_hit << "\n";
			}
		}

		offset_hit[g]->Write(TString::Format("offset_hit%d", g));
		offset_average_hit[g]->Write(TString::Format("offset_ave_hit%d", g));
		groups_hit->SetPoint(g, g, max_hit_offset);
		group_offset.push_back(max_hit_offset);
		last_offset = max_hit_offset;

		if (verbose_) {
			std::cout << "----------------------------------------\n"
				<< "group " << g << ", hit " << max_hit << ", offset " << max_hit_offset << "\n"
				<< "----------------------------------------\n";
		}


		// second children loop to save the search window
		group_window[g] = new TH1F(
			TString::Format("group_window%d", g), "window", 1000, -search_window_, search_window_
		);
		for (
			size_t vme_entry = group_size * g;
			vme_entry < group_size * (g + 1) && vme_entry < vme_times_.size();
			++vme_entry
		) {
			Long64_t vme_time = vme_times_[vme_entry];
			for (
				auto xia_iter = std::lower_bound(
					xia_times_.begin(),
					xia_times_.end(),
					vme_time - search_window_ + max_hit_offset
				);
				xia_iter != xia_times_.end();
				++xia_iter
			) {
				if (*xia_iter > vme_time + search_window_ + max_hit_offset) break;
				group_window[g]->Fill(*xia_iter - vme_time - max_hit_offset);
			}
		}
		group_window[g]->Write(TString::Format("group_window%d", g));


		// third children loop to fill the xia-timestamp:vme-timestamp calibration graph
		int result_point = 0;
		for (
			size_t vme_entry = group_size * g;
			vme_entry < group_size * (g + 1) && vme_entry < vme_times_.size();
			++vme_entry
		) {
			Long64_t vme_time = vme_times_[vme_entry];
			double xia_time;
			int local_hit = 0;
			for (
				auto xia_iter = std::lower_bound(
					xia_times_.begin(),
					xia_times_.end(),
					vme_time - search_window_ + group_offset[g]
				);
				xia_iter != xia_times_.end();
				++xia_iter
			) {
				if (*xia_iter > vme_time + search_window_ + group_offset[g]) break;
				++local_hit;
				xia_time = *xia_iter;
			}
			if (local_hit == 1) {
				result->SetPoint(result_point++, vme_time, xia_time);
			}
		}
	}
	groups_hit->Write("groups_hit");

	return result;
}


int Alignment::BuildResult(double *calibration_param) {
	// time window for all events
	TH1F *time_window =new TH1F(
		"ht", "total time window", 1000, -search_window_, search_window_
	);
	// result tree
	TTree *opt = new TTree("tree", "alignment result");
	// output data
	// output match VME trigger time recorded in XIA
	long long xia_time;
	// output vme trigger time recorded in VME
	long long vme_time;
	// setup branches
	opt->Branch("xia_time", &xia_time, "xt/L");
	opt->Branch("vme_time", &vme_time, "vt/L");


	AlignStatistics statistics(
		run_,
		xia_times_.size(),
		vme_times_.size(),
		calibration_param
	);

	// total entries in loop
	size_t entries = vme_times_.size();
	// 1/100 of total entries, for showing process
	size_t entry100 = entries / 100 + 1;
	// show begin
	printf("Building result   0%%");
	fflush(stdout);
	for (size_t vme_entry = 0; vme_entry < entries; ++vme_entry) {
		// show process
		if (vme_entry % entry100 == 0) {
			printf("\b\b\b\b%3ld%%", vme_entry / entry100);
			fflush(stdout);
		}

		// number of VME trigger in XIA match this VME event
		size_t match_count = 0;


		// calculate VME events time
		vme_time =
			calibration_param[0] +
			calibration_param[1] * vme_times_[vme_entry];

		// search for VME trigger in XIA
		for (
			auto xia_iter = std::lower_bound(
				xia_times_.begin(),
				xia_times_.end(),
				vme_time - search_window_
			);
			xia_iter != xia_times_.end();
			++xia_iter
		) {
			if (*xia_iter > vme_time + search_window_) break;
			++match_count;
			xia_time = *xia_iter;
		}
		if (match_count == 1) {
			++statistics.align_events;
			time_window->Fill(xia_time - vme_time);
		} else {
			++statistics.oversize_events;
			// set to -1.0 as placeholder
			xia_time = -1.0;
		}

		opt->Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	opt->Write();
	time_window->Write();

	statistics.Write();
	statistics.Print();

	return 0;
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

	// time offset, get from show_gdc_offset
	long long offset = run_ + 168;

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

	// second loop to read from first input file (this run)
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
	TString cache_file_name;
	cache_file_name.Form("%s%salign-%04d.root", kGenerateDataPath, kAlignDir, run_);
	TFile *cache_file = new TFile(cache_file_name, "recreate");

	TGraph *xia_vme_time = GroupAlignment();
	if (!xia_vme_time) {
		std::cerr << "Error: group alignment failed.\n";
		return -1;
	}

	// calibration
	TF1 *time_fit = new TF1("time_fit", "pol1", 0, 5e12);
	if (verbose_) {
		xia_vme_time->Fit(time_fit, "R+");
	} else {
		xia_vme_time->Fit(time_fit, "RQ+");
	}
	Double_t calibration_param[2];
	time_fit->GetParameters(calibration_param);
	std::cout << "Calibration p0 " << std::setprecision(15) << calibration_param[0]
		<< ", p1 " << std::setprecision(15) << calibration_param[1] << "\n";
	xia_vme_time->Write("calibration");

	BuildResult(calibration_param);

	cache_file->Close();


	// align gdc
	if (AlignGdc()) {
		return -1;
	}
	return 0;
}


}