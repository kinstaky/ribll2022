#include "include/alignment.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>

#include <TGraph.h>
#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>

#include "include/defs.h"

namespace ribll {

Alignment::Alignment(
	int run,
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
{
	verbose_ = false;
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
	Long64_t timestamp;
	ipt->SetBranchAddress("timestamp", &timestamp);

	// show process
	printf("reading xia events   0%%");
	fflush(stdout);
	Long64_t nentry100 = ipt->GetEntries() / 100;
	xia_times_.clear();
	for (Long64_t entry = 0; entry < ipt->GetEntries(); ++entry) {
		if (entry % nentry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / nentry100);
			fflush(stdout);
		}
		
		ipt->GetEntry(entry);
		xia_times_.push_back(timestamp);

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
	file_name.Form("%s%s%04d.root", kCrate4Path, kCrate4FileName, run_);
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
	Long64_t nentry100 = ipt->GetEntries() / 100;
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
		
			offset_hit[g]->AddPoint(offset, hit);
			Double_t average_hit = hit / ievent_count;
			if (max_hit < average_hit) {
				max_hit = average_hit;
				max_hit_offset = offset;
			}
			offset_average_hit[g]->AddPoint(offset, average_hit);


			if (verbose_ && (Long64_t)offset % 1'000'000'000 == 0) {
				std::cout << "  " << std::setprecision(10)
					<< offset / 1'000'000'000.0 << " s  "
					<< hit << "  " << average_hit << "\n";
			}
		}

		offset_hit[g]->Write(TString::Format("offset_hit%d", g));
		offset_average_hit[g]->Write(TString::Format("offset_ave_hit%d", g));
		groups_hit->AddPoint(g, max_hit_offset);
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
		for (
			size_t vme_entry = group_size * g;
			vme_entry < group_size * (g + 1) && vme_entry < vme_times_.size();
			++vme_entry
		) {
			Long64_t vme_time = vme_times_[vme_entry];
			Long64_t xia_time;
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
				result->AddPoint(vme_time, xia_time);
			}
		}
	}
	groups_hit->Write("groups_hit");

	return result;
}


int Alignment::BuildResult(Double_t *calibration_param) {
	// time window for all events
	TH1F *time_window = new TH1F("ht", "total time window", 1000, -search_window_, search_window_);
	// result tree
	TTree *opt = new TTree("tree", "alignment result");
	// output data
	Long64_t xia_time;
	Long64_t vme_time;
	// setup branches
	opt->Branch("xia_time", &xia_time, "xt/L");
	opt->Branch("vme_time", &vme_time, "vt/L");
	
	// file to record aligned timestamp of vme
	TString output_file_name;
	output_file_name.Form("%s%stimestamp-%04d.txt", kGenerateDataPath, kAlignDir, run_);
	std::ofstream fout(output_file_name.Data());
	if (!fout.good()) {
		std::cerr << "Error: open file " << output_file_name << " failed.\n";
		return -1;
	}

	// show process
	printf("Building result   0%%");
	fflush(stdout);
	size_t nentry100 = vme_times_.size() / 100;
	size_t align_events = 0;
	for (size_t vme_entry = 0; vme_entry < vme_times_.size(); ++vme_entry) {
		if (vme_entry % nentry100 == 0) {
			printf("\b\b\b\b%3ld%%", vme_entry / nentry100);
			fflush(stdout);
		}
		
		vme_time = calibration_param[0] + calibration_param[1] * vme_times_[vme_entry];
		size_t local_hit = 0;
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
			++local_hit;
			xia_time = *xia_iter;
		}
		if (local_hit == 1) {
			++align_events;
			time_window->Fill(xia_time - vme_time);
			opt->Fill();
		}

		fout << (local_hit == 1 ? xia_time : -1) << "\n";
	}
	printf("\b\b\b\b100%%\n");

	fout.close();
	opt->Write();
	time_window->Write();

	std::cout << "xia alignment rate " << align_events << " / " << xia_times_.size()
		<< " " << double(align_events) / double(xia_times_.size()) << "\n"
		<< "vme alignment rate " << align_events << " / " << vme_times_.size()
		<< " " << double(align_events) / double(vme_times_.size()) << "\n";

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
	return 0;
}

}