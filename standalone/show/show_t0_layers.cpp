#include <iostream>
#include <map>

#include <TFile.h>
#include <TString.h>

#include "include/event/dssd_event.h"

using namespace ribll;

struct LayerCount {
	unsigned int layer;
	unsigned int count;
};

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " run\n"
			<< "  run               Set run number.\n";
		return -1;
	}
	unsigned int run = atoi(argv[1]);

	// T0D1 merge file name
	TString d1_file_name;
	d1_file_name.Form(
		"%s%st0d1-merge-ta-%04u.root",
		kGenerateDataPath, kMergeDir, run
	);
	// T0D1 file
	TFile d1_file(d1_file_name, "read");
	// tree
	TTree *ipt = (TTree*)d1_file.Get("tree");
	if (!ipt) {
		std::cout << "Error: Get tree from "
			<< d1_file_name << " failed.\n";
		return -1;
	}
	// T0D2 merge file name
	TString d2_file_name;
	d2_file_name.Form(
		"%s%st0d2-merge-ta-%04u.root",
		kGenerateDataPath, kMergeDir, run
	);
	ipt->AddFriend("d2=tree", d2_file_name);
	// T0D3 merge file name
	TString d3_file_name;
	d3_file_name.Form(
		"%s%st0d3-merge-ta-%04u.root",
		kGenerateDataPath, kMergeDir, run
	);
	ipt->AddFriend("d3=tree", d3_file_name);
	// d1 event
	DssdMergeEvent d1_event;
	d1_event.SetupInput(ipt);
	// d2 event
	DssdMergeEvent d2_event;
	d2_event.SetupInput(ipt, "d2.");
	// d3 event
	DssdMergeEvent d3_event;
	d3_event.SetupInput(ipt, "d3.");

	LayerCount layer_counts[27];
	for (size_t i = 0; i < 27; ++i) {
		layer_counts[i].count = 0;
		layer_counts[i].layer = i;
	}

	// loop to get the merge event
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);
		unsigned int layer = d1_event.hit;
		layer += d2_event.hit * 3;
		layer += d3_event.hit * 9;
		layer_counts[layer].count += 1;
	}

	std::sort(
		layer_counts,
		layer_counts+27,
		[](const LayerCount &x, const LayerCount &y) {
			return x.count > y.count;
		}
	);

	for (const auto &[layer, count]: layer_counts) {
		unsigned int decimal_layer = layer % 3;
		decimal_layer += (layer / 3) % 3 * 10;
		decimal_layer += (layer / 9) * 100;
		std::cout << decimal_layer << "  " << count << "\n";
	}

	return 0;
}