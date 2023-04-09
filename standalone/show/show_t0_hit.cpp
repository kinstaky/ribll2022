#include <iostream>
#include <string>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include "include/defs.h"
#include "include/statistics/t0_hit_statistics.h"

using namespace ribll;

void PrintHit(
	unsigned int run,
	const std::string &tag
) {
	// t0d1 file name
	TString t0d1_file_name;
	t0d1_file_name.Form(
		"%s%st0d1-fundamental-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		tag.empty() ? "" : (tag + "-").c_str(),
		run
	);
	// file
	TFile ipf(t0d1_file_name, "read");
	// tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< t0d1_file_name << " failed.\n";
		return;
	}
	// t0d2 file name
	TString t0d2_file_name;
	t0d2_file_name.Form(
		"%s%st0d2-fundamental-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		tag.empty() ? "" : (tag + "-").c_str(),
		run
	);
	ipt->AddFriend("d2=tree", t0d2_file_name);
	// input data
	unsigned short d1_front_hit;
	unsigned short d2_front_hit;
	// setup branch
	ipt->SetBranchAddress("front_hit", &d1_front_hit);
	ipt->SetBranchAddress("d2.front_hit", &d2_front_hit);

	T0HitStatistics statistics(run, tag, ipt->GetEntries());

	// total entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Checking hit rate   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		ipt->GetEntry(entry);
		if (d1_front_hit == 0 || d2_front_hit == 0) continue;
		statistics.d1_single_hit += d1_front_hit == 1 ? 1 : 0;
		statistics.d1_multi_hit += d1_front_hit > 1 ? 1 : 0;
		statistics.d2_single_hit += d2_front_hit == 1 ? 1 : 0;
		statistics.d2_multi_hit += d2_front_hit > 1 ? 1 : 0;
		statistics.d1d2_single_hit += d1_front_hit > 1 && d2_front_hit > 1 ? 0 : 1;
		statistics.d1d2_multi_hit += d1_front_hit > 1 && d2_front_hit > 1 ? 1 : 0;
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// close file
	ipf.Close();

	statistics.Write();
	statistics.Print();

	return;
}


int main(int argc, char **argv) {
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " run\n"
			<< "  run               Set run number.\n";
		return -1;
	}

	unsigned int run = atoi(argv[1]);

	PrintHit(run, "ta");
	PrintHit(run, "");

	return 0;
}