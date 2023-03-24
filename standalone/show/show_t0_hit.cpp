#include <iostream>
#include <string>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include "include/defs.h"

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
	TFile *ipf = new TFile(t0d1_file_name, "read");
	// tree
	TTree *ipt = (TTree*)ipf->Get("tree");
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


	long long d1_single_hit = 0;
	long long d1_multi_hit = 0;
	long long d2_single_hit = 0;
	long long d2_multi_hit = 0;
	long long d1d2_single_hit = 0;
	long long d1d2_multi_hit = 0;

	long long entries = ipt->GetEntries();
	for (long long entry = 0; entry < entries; ++entry) {
		ipt->GetEntry(entry);
		if (d1_front_hit == 0 || d2_front_hit == 0) continue;
		d1_single_hit += d1_front_hit == 1 ? 1 : 0;
		d1_multi_hit += d1_front_hit > 1 ? 1 : 0;
		d2_single_hit += d2_front_hit == 1 ? 1 : 0;
		d2_multi_hit += d2_front_hit > 1 ? 1 : 0;
		d1d2_single_hit += d1_front_hit > 1 && d2_front_hit > 1 ? 0 : 1;
		d1d2_multi_hit += d1_front_hit > 1 && d2_front_hit > 1 ? 1 : 0;
	}

	std::cout << run << "," << tag << "," << entries
		<< "," << d1_single_hit << "," << d1_multi_hit
		<< "," << d2_single_hit << "," << d2_multi_hit
		<< "," << d1d2_single_hit << "," << d1d2_multi_hit
		<< "," << double(d1_single_hit) / double(d1_multi_hit)
		<< "," << double(d2_single_hit) / double(d2_multi_hit)
		<< "," << double(d1d2_single_hit) / double(d1d2_multi_hit)
		<< "\n";
	return;
}


int main(int argc, char **argv) {
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " run\n"
			<< "  run               Set run number.\n";
		return -1;
	}

	unsigned int run = atoi(argv[1]);

	PrintHit(run, "");
	PrintHit(run, "ta");


	return 0;
}