#include "include/time_reference.h"

#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>

namespace ribll {

int TimeReference(unsigned int run, const std::string &tag) {
	// input tof file
	TString tof_file_name;
	tof_file_name.Form(
		"%s%stof-fundamental-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// ppac file
	TFile ipf(tof_file_name, "read");
	// ppac tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< tof_file_name << " failed.\n";
		return -1;
	}
	// input xppac file
	// TString ppac_file_name;
	// ppac_file_name.Form(
	// 	"%s%sxppac-fundamental-%s%04u.root",
	// 	kGenerateDataPath,
	// 	kFundamentalDir,
	// 	tag.empty() ? "" : (tag+"-").c_str(),
	// 	run
	// );
	// // add friend
	// ipt->AddFriend("ppac=tree", ppac_file_name);
	// input tof event
	TofFundamentalEvent tof_event;
	// // input xppac event
	// PpacFundamentalEvent ppac_event;
	// setup input branches
	tof_event.SetupInput(ipt);
	// ppac_event.SetupInput(ipt, "ppac.");

	// output reference file name
	TString ref_file_name;
	ref_file_name.Form(
		"%s%sreftime-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// output reference file
	TFile ref_file(ref_file_name, "recreate");
	// output reference tree
	TTree ref_tree("tree", "reference time");
	// output reference time
	double ref_time;
	// setup output tree branch
	ref_tree.Branch("time", &ref_time, "t/D");

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Calculating reference time   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// read event
		ipt->GetEntry(entry);
		// initialize
		ref_time = -1e5;
		// t1 and t2
		if (
			tof_event.time[0] > -9e4
			&& tof_event.time[1] > -9e4
			&& tof_event.time[0]-tof_event.time[1] > -60
			&& tof_event.time[0]-tof_event.time[1] < -57
		) {
			ref_time = tof_event.time[1];
		}

		// if (
		// 	(ppac_event.flag & 0x4200) == 0x4200
		// 	&& ppac_event.anode[1] - ppac_event.anode[2] > -16
		// 	&& ppac_event.anode[1] - ppac_event.anode[2] < -13
		// ) {
		// 	// a1 and a2
		// 	ref_time = ppac_event.anode[2];
		// } else if (
		// 	(ppac_event.flag & 0x4010) == 0x4010
		// 	&& ppac_event.anode[0] - ppac_event.anode[2] > -18
		// 	&& ppac_event.anode[0] - ppac_event.anode[2] < -13
		// ) {
		// 	// a0 and a2
		// 	ref_time = ppac_event.anode[2];
		// } else if (
		// 	(ppac_event.flag & 0x210) == 0x210
		// 	&& ppac_event.anode[0] - ppac_event.anode[1] > -2
		// 	&& ppac_event.anode[0] - ppac_event.anode[1] < 1.5
		// ) {
		// 	// a0 and a1
		// 	ref_time = ppac_event.anode[1] + 14.2;
		// } else if (
		// 	(ppac_event.flag & 0x4000) == 0x4000
		// 	&& tof_event.time[1] > -9e4
		// 	&& tof_event.time[1] - ppac_event.anode[2] > -2
		// 	&& tof_event.time[1] - ppac_event.anode[2] < 5
		// ) {
		// 	// t1 and a2
		// 	ref_time = ppac_event.anode[2];
		// } else if (
		// 	(ppac_event.flag & 0x200) == 0x200
		// 	&& tof_event.time[1] > -9e4
		// 	&& tof_event.time[1] - ppac_event.anode[1] > 13
		// 	&& tof_event.time[1] - ppac_event.anode[1] < 17
		// ) {
		// 	// t1 and a1
		// 	ref_time = ppac_event.anode[1] + 14.2;
		// } else if (
		// 	(ppac_event.flag & 0x10) == 0x10
		// 	&& tof_event.time[1] > -9e4
		// 	&& tof_event.time[1] - ppac_event.anode[0] > 13
		// 	&& tof_event.time[1] - ppac_event.anode[0] < 17
		// ) {
		// 	// t1 and a0
		// 	ref_time = ppac_event.anode[0] + 13.9;
		// } else if (
		// 	(ppac_event.flag & 0x4000) == 0x4000
		// 	&& tof_event.time[0] > -9e4
		// 	&& tof_event.time[0] - ppac_event.anode[2] > -60
		// 	&& tof_event.time[0] - ppac_event.anode[2] < -55
		// ) {
		// 	// t0 and a2
		// 	ref_time = ppac_event.anode[2];
		// } else if (
		// 	(ppac_event.flag & 0x200) == 0x200
		// 	&& tof_event.time[0] > -9e4
		// 	&& tof_event.time[0] - ppac_event.anode[1] > -46
		// 	&& tof_event.time[0] - ppac_event.anode[1] < -41
		// ) {
		// 	// t0 and a1
		// 	ref_time = ppac_event.anode[1] + 14.2;
		// } else if (
		// 	(ppac_event.flag & 0x10) == 0x10
		// 	&& tof_event.time[0] > -9e4
		// 	&& tof_event.time[0] - ppac_event.anode[0] > -46
		// 	&& tof_event.time[0] - ppac_event.anode[0] < -41
		// ) {
		// 	// t0 and a0
		// 	ref_time = ppac_event.anode[0] + 13.9;
		// }
		ref_tree.Fill();
	}
	printf("\b\b\b\b100%%\n");

	// save tree and close files
	ref_tree.Write();
	ref_file.Close();
	ipf.Close();

	return 0;
}

}		// time reference