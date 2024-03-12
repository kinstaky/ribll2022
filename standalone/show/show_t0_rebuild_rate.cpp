/*
 * 这个程序计算 T0 中重建出 4He 和 10Be 的比例，用于衡量重建算法（merge，track）的好坏。
 */

#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

#include "include/event/particle_event.h"

using namespace ribll;

int main(int argc, char **argv) {
	int run = 0;
	if (argc >= 2) {
		run = atoi(argv[1]);
	}

	// input T0 particle file name
	TString t0_file_name = TString::Format(
		"%s%st0-particle-sim-ta-%04d.root",
		kGenerateDataPath, kParticleDir, run
	);
	// input T0 particle file
	TFile ipf(t0_file_name, "read");
	// input T0 particle tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< t0_file_name << " failed.\n";
		return -1;
	}
	// input particle event
	bool has_4he, has_10be;
	// setup input branches
	ipt->SetBranchAddress("has_4he", &has_4he);
	ipt->SetBranchAddress("has_10be", &has_10be);

	// input T0 detect file name
	TString t0_detect_file_name = TString::Format(
		"%s%st0-detect-%04d.root",
		kGenerateDataPath, kSimulateDir, run
	);
	// add friend
	ipt->AddFriend("detect=tree", t0_detect_file_name);
	// T0 detect DSSD flags
	int dssd_flag[3];
	// T0 detect T0 flag
	int t0_flag;
	// energy cut by threshold times
	int cut_time;
	ipt->SetBranchAddress("dssd_flag", dssd_flag);
	ipt->SetBranchAddress("t0_flag", &t0_flag);
	ipt->SetBranchAddress("cut_time", &cut_time);

	// statistics
	// total count of different T0 flags
	// 0~5: t0_flag 0~5, 6: one cut, 7: total
	int total_count[8];

	// rebuild count of different T0 flags
	// 0~5: t0_flag 0~5, 6: one cut, 7: total
	int rebuild_count[8];

	// initialize
	for (int i = 0; i < 8; ++i) {
		total_count[i] = 0;
		rebuild_count[i] = 0;
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Calculating rebuild rate   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get data
		ipt->GetEntry(entry);

		if (cut_time == 1) {
			++total_count[6];
		} else {
			++total_count[t0_flag];
		}
		++total_count[7];
		if (has_10be && has_4he) {
			if (cut_time == 1) {
				++rebuild_count[6];
			} else {
				++rebuild_count[t0_flag];
			}
			++rebuild_count[7];

		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// close files
	ipf.Close();

	// show statistics
	std::cout << "Impossible " << rebuild_count[0] << " / " << total_count[0]
		<< " " << double(rebuild_count[0])/double(total_count[0]) << "\n"
		<< "Only 2 " << rebuild_count[1] << " / " << total_count[1]
		<< " " << double(rebuild_count[1])/double(total_count[1]) << "\n"
		<< "Include 2a-2a " << rebuild_count[2] << " / " << total_count[2]
		<< " " << double(rebuild_count[2])/double(total_count[2]) << "\n"
		<< "Include 1-2 " << rebuild_count[3] << " / " << total_count[3]
		<< " " << double(rebuild_count[3])/double(total_count[3]) << "\n"
		<< "Include 1-2a " << rebuild_count[4] << " / " << total_count[4]
		<< " " << double(rebuild_count[4])/double(total_count[4]) << "\n"
		<< "Include 1-1 " << rebuild_count[5] << " / " << total_count[5]
		<< " " << double(rebuild_count[5])/double(total_count[5]) << "\n"
		<< "One cut " << rebuild_count[6] << " / " << total_count[6]
		<< " " << double(rebuild_count[6])/double(total_count[6]) << "\n"
		<< "Total " << rebuild_count[7] << " / " << total_count[7]
		<< " " << double(rebuild_count[7])/double(total_count[7]) << "\n";
	return 0;
}