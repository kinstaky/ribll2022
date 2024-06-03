#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TH2F.h>

#include "include/defs.h"

using namespace ribll;

int main() {
	// total number of run
	int total = 84;
	// number of run for pixel regard as hole
	int hole[32][32];
	for (int i = 0; i < 32; ++i) {
		for (int j = 0; j < 32; ++j) {
			hole[i][j] = 0;
		}
	}

	// input file name
	TString input_file_name = TString::Format(
		"%s%shole-flag.root", kGenerateDataPath, kHoleDir
	);
	// input file
	TFile ipf(input_file_name, "read");
	// get histogram
	TH2F *hhsr = (TH2F*)ipf.Get("hhsr");
	if (!hhsr) {
		std::cerr << "Error: Get histogram from "
			<< input_file_name << " failed.\n";
		return -1;
	}

	for (int fs = 0; fs < 32; ++fs) {
		for (int bs = 0; bs < 32; ++bs) {
			int start_run = hhsr->GetBinContent(fs+1, bs+1);
			if (start_run == 0) continue;
			for (int run = 618; run <= 746; ++run) {
				if (run == 628) continue;
				if (run > 652 && run < 675) continue;
				if (run > 716 && run < 739) continue;
				if (run >= start_run) ++hole[fs][bs];
			}
		}
	}

	// output file name
	TString output_file_name = TString::Format(
		"%s%saverage-hole.root",
		kGenerateDataPath, kHoleDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output histogram, average hole possibility
	TH2F hahp("hahp", "possibility to become hole", 32, 0, 32, 32, 0, 32);
	for (int i = 0; i < 32; ++i) {
		for (int j = 0; j < 32; ++j) {
			hahp.SetBinContent(i+1, j+1, double(hole[i][j])/double(total));
		}
	}
	hahp.Write();
	opf.Close();


	return 0;
}