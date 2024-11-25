#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TGraph.h>

#include "include/event/threebody_info_event.h"
#include "include/event/channel_v2_event.h"

using namespace ribll;

// straight_cut parameters
constexpr double sa12 = 0.47;
constexpr double sb12 = -0.042;
constexpr double sa23 = 0.59;
constexpr double sb23 = -0.042;


inline double Straight(double de, double e, double a, double b) {
	return sqrt(de*e + a*de*de) + b*e;
}

int main(int argc, char **argv) {
	if (argc != 2) {
		std::cout << "Usage: " << argv[0] << " run\n"
			"  run            run number\n";
		return -1;
	}
	// run number
	int run = atoi(argv[1]);
	if (run < 1 || run > 2) {
		std::cerr << "Error: Invalid run " << run << "\n";
		return -1;
	}

	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody-sim-%04d.root", kGenerateDataPath, kInformationDir, run
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
	// input event
	ThreeBodyInfoEvent event;
	int straight_cut;
	// setup input branches
	event.SetupInput(ipt);

	// output file name
	TString output_file_name = TString::Format(
		"%s%sefficiency-%04d.root",
		kGenerateDataPath,
		kSimulateDir,
		run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// T0 position efficiency
	TGraph graph_t0_position_efficiency[3];
	// T0 not hole efficiency
	TGraph graph_t0_not_hole_efficiency[3];
	// T0 not bind efficiency
	TGraph graph_t0_not_bind_efficiency[3];
	// T0 in straight PID cut efficiency
	TGraph graph_t0_straight_cut_efficiency[3];
	// T0 efficiency graph
	TGraph graph_t0_efficiency[3];
	// TAF efficiency graph
	TGraph graph_taf_efficiency[3];
	// XIA PPAC efficiency
	TGraph graph_xppac_efficiency[3];
	// target efficiency
	TGraph graph_target_efficiency[3];
	// total efficiency
	TGraph graph_efficiency[3];

	int t0_position_count[3][100];
	int t0_not_hole_count[3][100];
	int t0_not_bind_count[3][100];
	int t0_straight_cut_count[3][100];
	int t0_valid_count[3][100];
	int taf_valid_count[3][100];
	int xppac_valid_count[3][100];
	int target_valid_count[3][100];
	int valid_count[3][100];
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 100; ++j) {
			t0_position_count[i][j] = 0;
			t0_not_hole_count[i][j] = 0;
			t0_not_bind_count[i][j] = 0;
			t0_straight_cut_count[i][j] = 0;
			t0_valid_count[i][j] = 0;
			taf_valid_count[i][j] = 0;
			xppac_valid_count[i][j] = 0;
			target_valid_count[i][j] = 0;
			valid_count[i][j] = 0;
		}
	}

	long long entries = ipt->GetEntries();
	// 1/100 of total entries, for showing process
	long long entry100 = entries / 100;
	// show start
	printf("Getting efficiency   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get data
		ipt->GetEntry(entry);

		// index
		int i1 = entry / (entries / 3);
		int i2 = (entry % (entries / 3)) / (entries / 300);

		// calculate straight_cut flag
		straight_cut = 0;
		if (event.layer[0] == 1) {
			double ef = Straight(
				event.be_channel[0], event.be_channel[1], sa12, sb12
			);
			if (ef > 20200 && ef < 21200) straight_cut |= 1;
		} else {
			double ef = Straight(
				event.be_channel[1], event.be_channel[2], sa23, sb23
			);
			if (ef > 24000 && ef < 25000) straight_cut |= 1;
		}
		if (event.layer[1] == 1) {
			double ef = Straight(
				event.he_channel[0], event.he_channel[1], sa12, sb12
			);
			if (ef > 6100 && ef < 6600) straight_cut |= 2;
		} else if (event.layer[1] == 2) {
			double ef = Straight(
				event.he_channel[1], event.he_channel[2], sa23, sb23
			);
			if (ef > 7000 && ef < 8000) straight_cut |= 2;
		} else {
			straight_cut |= 2;
		}

		// T0 position valid
		if (event.taf_flag == 0 || event.taf_flag == 5) {
			++t0_position_count[i1][i2];
		}
		if (!event.hole[0] && !event.hole[1]) {
			++t0_not_hole_count[i1][i2];
		}
		if (event.bind == 0) {
			++t0_not_bind_count[i1][i2];
		}
		if (straight_cut == 3) {
			++t0_straight_cut_count[i1][i2];
		}

		// T0 valid
		if (
			(event.taf_flag == 0 || event.taf_flag == 5)
			&& !(event.hole[0] || event.hole[1])
			&& event.bind == 0
			&& straight_cut == 3
		) {
			++t0_valid_count[i1][i2];
		}
		// TAF valid
		if (event.taf_flag == 0 || event.taf_flag == 4) {
			++taf_valid_count[i1][i2];
		}
		// PPAC valid
		if (event.xppac_xflag != 0 && event.xppac_yflag != 0) {
			++xppac_valid_count[i1][i2];
		}
		// target valid
		if (event.target_flag == 1) {
			++target_valid_count[i1][i2];
		}
		// All valid
		if (
			!(event.hole[0] || event.hole[1])
			&& event.bind == 0
			&& straight_cut == 3
			&& event.taf_flag == 0
			&& (event.xppac_xflag != 0 && event.xppac_yflag != 0)
			&& event.target_flag == 1
		) {
			++valid_count[i1][i2];
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// close file
	ipf.Close();


	// fill channel-V2 events
	TString v2_file_name = TString::Format(
		"%s%sC14-10Be-4He-2H-v2-sim.root",
		kGenerateDataPath,
		kChannelDir
	);
	// input file
	TFile v2_ipf(v2_file_name, "read");
	// input tree
	TTree *v2_ipt = (TTree*)v2_ipf.Get("tree");
	if (!v2_ipt) {
		std::cerr << "Error: Get tree from "
			<< v2_file_name << " failed.\n";
		return -1;
	}
	// input event
	ChannelV2Event v2_event;
	// setup input branches
	v2_event.SetupInput(v2_ipt);
	// show start
	printf("Getting v2 efficiency   0%%");
	fflush(stdout);
	// total number of entries
	long long v2_entries = v2_ipt->GetEntries();
	entry100 = v2_entries / 100 + 1;
	for (long long entry = 0; entry < v2_entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get data
		v2_ipt->GetEntry(entry);

		// index
		int i1 = v2_event.entry / (entries / 3);
		int i2 = (v2_event.entry % (entries / 3)) / (entries / 300);

		// TAF valid
		if (!v2_event.tafcsi_valid && v2_event.tafd_front_strip >= 13) {
			++taf_valid_count[i1][i2];
		}

		int target_flag = 0;
		if (pow(v2_event.tx-2.5, 2.0) + pow(v2_event.ty+2.0, 2.0) < 196.0) {
			target_flag = 1;
		}
		// All valid
		if (
			v2_event.valid == 0
			&& !v2_event.tafcsi_valid
			&& v2_event.tafd_front_strip >= 13
			&& v2_event.ppac_xnum >= 1
			&& v2_event.ppac_ynum >= 1
			&& target_flag == 1
		) {
			++valid_count[i1][i2];
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	// close file
	v2_ipf.Close();

	const double energy_base[3] = {12.02, 15.39, 18.20};
	// fill to graph
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 100; ++j) {
			double energy = energy_base[i] + 0.2 * j;
			graph_t0_position_efficiency[i].AddPoint(
				energy, t0_position_count[i][j] / double(entries / 300)
			);
			graph_t0_not_hole_efficiency[i].AddPoint(
				energy, t0_not_hole_count[i][j] / double(entries / 300)
			);
			graph_t0_not_bind_efficiency[i].AddPoint(
				energy, t0_not_bind_count[i][j] / double(entries / 300)
			);
			graph_t0_straight_cut_efficiency[i].AddPoint(
				energy, t0_straight_cut_count[i][j] / double(entries / 300)
			);
			graph_t0_efficiency[i].AddPoint(
				energy, t0_valid_count[i][j] / double(entries / 300)
			);
			graph_taf_efficiency[i].AddPoint(
				energy, taf_valid_count[i][j] / double(entries / 300)
			);
			graph_xppac_efficiency[i].AddPoint(
				energy, xppac_valid_count[i][j] / double(entries / 300)
			);
			graph_target_efficiency[i].AddPoint(
				energy, target_valid_count[i][j] / double(entries / 300)
			);
		}
		graph_efficiency[i].AddPoint(
			energy_base[i],
			valid_count[i][0] / double(entries / 300)
		);
		for (int j = 1; j < 99; ++j) {
			graph_efficiency[i].AddPoint(
				energy_base[i] + 0.2 * j,
				(valid_count[i][j-1] + valid_count[i][j] + valid_count[i][j+1])
					/ double(entries / 100)
			);
		}
		graph_efficiency[i].AddPoint(
			energy_base[i] + 0.2 * 99,
			valid_count[i][99] / double(entries / 300)
		);
	}

	// save graphs
	opf.cd();
	for (int i = 0; i < 3; ++i) {
		graph_t0_position_efficiency[i].Write(TString::Format("t0pg%d", i));
	}
	for (int i = 0; i < 3; ++i) {
		graph_t0_not_hole_efficiency[i].Write(TString::Format("t0hg%d", i));
	}
	for (int i = 0; i < 3; ++i) {
		graph_t0_not_bind_efficiency[i].Write(TString::Format("t0bg%d", i));
	}
	for (int i = 0; i < 3; ++i) {
		graph_t0_straight_cut_efficiency[i].Write(
			TString::Format("t0sg%d", i)
		);
	}
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
		graph_target_efficiency[i].Write(TString::Format("targetg%d", i));
	}
	for (int i = 0; i < 3; ++i) {
		graph_efficiency[i].Write(TString::Format("g%d", i));
	}
	// close file
	opf.Close();
	return 0;
}
