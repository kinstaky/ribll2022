#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TH2F.h>

#include "include/defs.h"

using namespace ribll;

int main() {
	// total number of run
	int total = 0;
	// number of run for pixel regard as hole
	int hole[32][32];
	for (int i = 0; i < 32; ++i) {
		for (int j = 0; j < 32; ++j) {
			hole[i][j] = 0;
		}
	}

	// loop to check hole
	for (int run = 618; run < 746; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		if (run > 716 && run < 739) continue;
		// input file name
		TString input_file_name = TString::Format(
			"%s%sshow-t0d2-pixel-%04d.root",
			kGenerateDataPath, kShowDir, run
		);
		// input file
		TFile ipf(input_file_name, "read");
		// get histogram
		TH2F *hbr = (TH2F*)ipf.Get("hbr");
		if (!hbr) {
			std::cerr << "Error: Get histogram from "
				<< input_file_name << " failed.\n";
			return -1;
		}
		for (int i = 0; i < 32; ++i) {
			for (int j = 0; j < 32; ++j) {
				if (hbr->GetBinContent(i+1, j+1) > 0.08) {
					++hole[i][j];
				}
			}
		}
		ipf.Close();
		++total;
	}

	// output file name
	TString output_file_name = TString::Format(
		"%s%st0d2-average-hole.root",
		kGenerateDataPath, kShowDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output histogram, average hole possibility
	TH2F hahp("hahp", "average hole possibility", 32, 0, 32, 32, 0, 32);
	for (int i = 0; i < 32; ++i) {
		for (int j = 0; j < 32; ++j) {
			hahp.SetBinContent(i+1, j+1, double(hole[i][j])/double(total));
		}
	}
	hahp.Write();
	opf.Close();


	return 0;
}