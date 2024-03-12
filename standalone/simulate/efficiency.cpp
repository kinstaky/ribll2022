#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TGraph.h>

#include "include/event/channel_event.h"

using namespace ribll;

int main() {
	// input channel file name
	TString input_file_name = TString::Format(
		"%s%sC14-10Be-4He-2H-sim-0001.root",
		kGenerateDataPath,
		kChannelDir
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}

	// input XIA PPAC file
	TString xppac_file_name = TString::Format(
		"%s%sxppac-particle-sim-ta-0001.root",
		kGenerateDataPath,
		kParticleDir
	);
	// add friend
	ipt->AddFriend("xppac=tree", xppac_file_name);

	// input channel event
	ChannelEvent channel;
	// XIA PPAC flag
	unsigned short xppac_xflag, xppac_yflag;
	// setup input branches
	channel.SetupInput(ipt);
	ipt->SetBranchAddress("xppac.xflag", &xppac_xflag);
	ipt->SetBranchAddress("xppac.yflag", &xppac_yflag);


	// output file name
	TString output_file_name = TString::Format(
		"%s%sefficiency-0001.root",
		kGenerateDataPath,
		kSimulateDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// T0 efficiency graph
	TGraph graph_t0_efficiency[3];
	// TAF efficiency graph
	TGraph graph_taf_efficiency[3];
	// XIA PPAC efficiency
	TGraph graph_xppac_efficiency[3];
	// total efficiency
	TGraph graph_efficiency[3];

	int t0_valid_count[3][100];
	int taf_valid_count[3][100];
	int xppac_valid_count[3][100];
	int valid_count[3][100];
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 100; ++j) {
			t0_valid_count[i][j] = 0;
			taf_valid_count[i][j] = 0;
			xppac_valid_count[i][j] = 0;
			valid_count[i][j] = 0;
		}
	}

	if (ipt->GetEntries() != 3'000'000) {
		std::cerr << "Error: entries != 3,000,000\n";
		return -1;
	}
	// 1/100 of total entries, for showing process
	long long entry100 = 30'000;
	// show start
	printf("Getting efficiency   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < 3'000'000; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get data
		ipt->GetEntry(entry);

		// index
		int i1 = entry / 1'000'000;
		int i2 = (entry % 1'000'000) / 10'000;
		if (xppac_xflag != 0 && xppac_yflag != 0) {
			++xppac_valid_count[i1][i2];
		}
		if (channel.taf_index >= 0 && channel.taf_index <=5) {
			++t0_valid_count[i1][i2];
			++taf_valid_count[i1][i2];
			if (xppac_xflag != 0 && xppac_yflag != 0) {
				++valid_count[i1][i2];
			}
		} else if (channel.taf_index == -2) {
			++taf_valid_count[i1][i2];
		} else if (channel.taf_index == -4) {
			++t0_valid_count[i1][i2];
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// for (int i = 0; i < 100; ++i) {
	// 	std::cout << t0_valid_count[0][i] << ", "
	// 		<< taf_valid_count[0][i] << ", "
	// 		<< xppac_valid_count[0][i] << ", "
	// 		<< valid_count[0][i] << "\n";
	// }

	const double energy_base[3] = {12.02, 15.39, 18.20};
	// fill to graph
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 100; ++j) {
			double energy = energy_base[i] + 0.2 * j;
			graph_t0_efficiency[i].AddPoint(
				energy, t0_valid_count[i][j] / 10000.0
			);
			graph_taf_efficiency[i].AddPoint(
				energy, taf_valid_count[i][j] / 10000.0
			);
			graph_xppac_efficiency[i].AddPoint(
				energy, xppac_valid_count[i][j] / 10000.0
			);
			graph_efficiency[i].AddPoint(
				energy, valid_count[i][j] / 10000.0
			);
		}
	}

	// save graphs
	for (int i = 0; i < 3; ++i) {
		graph_t0_efficiency[i].Write(TString::Format("t0g%d", i));
	}
	for (int i = 0; i < 3; ++i) {
		graph_taf_efficiency[i].Write(TString::Format("tafg%d", i));
	}
	for (int i = 0; i < 3; ++i) {
		graph_xppac_efficiency[i].Write(TString::Format("xppacg%d", i));
	}
	for (int i = 0; i < 3; ++i) {
		graph_efficiency[i].Write(TString::Format("g%d", i));
	}
	// close file
	opf.Close();
	ipf.Close();
	return 0;
}