#include <iostream>
#include <vector>
#include <fstream>

#include <TFile.h>
#include <TString.h>
#include <TH2F.h>

#include "include/defs.h"

using namespace ribll;


int main() {
	int run_index = 0;
	int run_number[100];
	bool hole[32][32][100];
	int start_run[32][32];
	for (int fs = 0; fs < 32; ++fs) {
		for (int bs = 0; bs < 32; ++bs) {
			start_run[fs][bs] = 2000;
		}
	}
	for (int run = 618; run <= 746; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		if (run > 716 && run < 739) continue; 
		// record run number
		run_number[run_index] = run;

		// input file name
		TString input_file_name = TString::Format(
			"%s%sshow-t0d2-pixel-%04d.root",
			kGenerateDataPath, kShowDir, run
		);
		// input file
		TFile ipf(input_file_name, "read");
		// resolution histogram
		TH2F *hbr = (TH2F*)ipf.Get("hbr");
		if (!hbr) {
			std::cerr << "Error: Get pixel resolution from "
				<< input_file_name << " failed.\n";
			return -1;
		}
		for (int fs = 0; fs < 32; ++fs) {
			for (int bs = 0; bs < 32; ++bs) {
				double res = hbr->GetBinContent(fs+1, bs+1);
				hole[fs][bs][run_index] = res > 0.08;
			}
		}
		// close files
		ipf.Close();
		++run_index;
	}
	std::cout << "Total " << run_index << "\n";

	// output file name
	TString output_file_name = TString::Format(
		"%s%shole-flag.root",
		kGenerateDataPath, kHoleDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// histogram of start run
	TH2F hist_hole_start_run("hhsr", "start run of hole", 32, 0, 32, 32, 0, 32);

	for (int fs = 0; fs < 32; ++fs) {
		for (int bs = 0; bs < 32; ++bs) {
			for (int i = 0; i < run_index; ++i) {
				if (run_number[i] == 631 || run_number[i] == 695) continue;
				if (hole[fs][bs][i]) {
					start_run[fs][bs] = run_number[i];
					break;
				}
			}
		}
	}

	// manual correct for 8% resolution
	start_run[9][19] = 690;
	start_run[9][20] = 741;
	start_run[11][14] = 703;
	start_run[11][15] = 629;
	start_run[12][14] = 651;

	// manual correct for 5% resolution
	// start_run[9][13] = 0;
	// start_run[9][19] = 649;
	// start_run[9][20] = 683;
	// start_run[10][13] = 0;
	// start_run[10][15] = 623;
	// start_run[11][13] = 0;
	// start_run[11][14] = 0;
	// start_run[11][15] = 618;
	// start_run[11][22] = 627;
	// start_run[12][13] = 0;
	// start_run[12][14] = 618;
	// start_run[12][22] = 648;
	// start_run[13][13] = 0;
	// start_run[13][14] = 630;
	// start_run[13][22] = 644;
	// start_run[14][13] = 0;
	// start_run[15][13] = 0;
	// start_run[15][22] = 675;
	// start_run[17][14] = 0;
	// start_run[17][20] = 638;
	// start_run[18][13] = 0;
	// start_run[18][15] = 0;
	// start_run[18][21] = 714;
	// start_run[19][13] = 0;


	// check
	for (int fs = 0; fs < 32; ++fs) {
		for (int bs = 0; bs < 32; ++bs) {
			if (start_run[fs][bs] == 2000) continue;
			for (int i = 0; i < run_index; ++i) {
				if (run_number[i] <= start_run[fs][bs]) continue;
				if (!hole[fs][bs][i]) {
					std::cout << fs << ", " << bs << ", "
						<< run_number[i] << "\n";
					break;
				}
			}
		}
	}

	for (int fs = 0; fs < 32; ++fs) {
		for (int bs = 0; bs < 32; ++bs) {
			if (start_run[fs][bs] == 2000) continue;
			hist_hole_start_run.SetBinContent(
				fs+1, bs+1, start_run[fs][bs]
			);
		}
	}

	// output csv file name
	TString csv_file_name = TString::Format(
		"%s%shole-flag.csv", kGenerateDataPath, kHoleDir
	);
	// output stream
	std::ofstream fout(csv_file_name.Data());
	fout << "fs, bs, run, hole\n";
	for (int fs = 0; fs < 32; ++fs) {
		for (int bs = 0; bs < 32; ++bs) {
			if (start_run[fs][bs] == 2000) continue;
			for (int i = 0; i < run_index; ++i) {
				fout << fs << ", " << bs << ", "
					<< run_number[i] << ", " << hole[fs][bs][i] << "\n";
			}
		}
	}
	// close output stream
	fout.close();

	// save histograms
	hist_hole_start_run.Write();
	// close files
	opf.Close();

	return 0;
}