#include <iostream>
#include <fstream>
#include <string>

#include <TFile.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TString.h>

#include "include/defs.h"

using namespace ribll;

int main(int argc, char **argv) {
	if (argc != 4) {
		std::cerr << "Error: Too few arguments\n";
		std::cout << "Usage: " << argv[0] << " start_run end_run name\n";
		return -1;
	}
	int start_run = atoi(argv[1]);
	int end_run = atoi(argv[2]);
	std::string name(argv[3]);

	if (end_run < start_run) {
		std::cerr << "Error: Invalid end_run\n";
		return -1;
	}
	bool vppac = false;
	if (name == "vppac") {
		vppac = true;
	} else if (name == "xppac") {
		vppac = false;
	} else {
		std::cerr << "Error: Invalid name\n";
		return -1;
	}

	// output file name
	TString output_file_name = TString::Format(
		"%s%sshow-%s-change.root", kGenerateDataPath, kShowDir, name.c_str()
	);
	// output file
	TFile opf(output_file_name, "recreate");
	TGraph gsum[18];
	TMultiGraph mgsum[6];
	TGraph gdiff[6];
	TMultiGraph mgdiff[2];

	for (int run = start_run; run <= end_run; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		// sum range file name
		TString sum_file_name = TString::Format(
			"%s%s%s-range-%04d.txt",
			kGenerateDataPath, kTimeDir, name.c_str(), run
		);
		// input file stream
		std::ifstream fin(sum_file_name.Data());
		if (!fin.good()) {
			std::cerr << "Error: Open file " << sum_file_name << " failed.\n";
			continue;
		}
		// range data
		double range_left, range_right, range_center;
		// read data
		for (int i = 0; i < 18; ++i) {
			if (vppac) {
				fin >> range_left >> range_right >> range_center;
			} else {
				fin >> range_left >> range_right;
				range_center = (range_left + range_right) / 2.0;
			}
			gsum[i].AddPoint(run, range_center);
		}
		for (int i = 0; i < 6; ++i) {
			int offset = (i / 2) * 6 + (i % 2);
			mgsum[i].Add(gsum+offset, "APL");
			mgsum[i].Add(gsum+offset+2, "PL");
			mgsum[i].Add(gsum+offset+4, "PL");
		}
		fin.close();

		// normalize file name
		TString normalize_file_name = TString::Format(
			"%s%sppac-normalize-%04d.txt",
			kGenerateDataPath, kNormalizeDir, run
		);
		// input file
		std::ifstream norm_fin(normalize_file_name.Data());
		if (!norm_fin.good()) {
			std::cerr << "Error: Open file "
				<< normalize_file_name << " failed.\n";
			continue;
		}
		// offset data
		double offset;
		for (int i = 0; i < 6; ++i) {
			norm_fin >> offset;
			gdiff[i].AddPoint(run, offset);
		}
		norm_fin.close();
		for (int i = 0; i < 2; ++i) {
			mgdiff[i].Add(gdiff+i*3, "APL");
			mgdiff[i].Add(gdiff+i*3+1, "PL");
			mgdiff[i].Add(gdiff+i*3+2, "PL");
		}
	}

	// save graph
	for (int i = 0; i < 18; ++i) {
		gsum[i].Write(TString::Format(
			"gs%c%da%d", "xy"[i%2], i/6, (i/2)%3
		));
	}
	for (int i = 0; i < 6; ++i) {
		int offset = (i / 2) * 6 + (i % 2);
		gsum[offset+2].SetLineColor(kBlue);
		gsum[offset+4].SetLineColor(kRed);
		mgsum[i].Write(TString::Format(
			"mgs%c%d", "xy"[i%2], (i/2)%3
		));
	}
	for (int i = 0; i < 6; ++i) {
		gdiff[i].Write(TString::Format(
			"gd%c%d", "xy"[i/3], i%3
		));
	}
	for (int i = 0; i < 2; ++i) {
		gdiff[i*3+1].SetLineColor(kBlue);
		gdiff[i*3+2].SetLineColor(kRed);
		mgdiff[i].Write(TString::Format("mgd%c", "xy"[i%2]));
	}
	// close file
	opf.Close();
	return 0;
}