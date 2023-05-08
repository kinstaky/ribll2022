#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TString.h>

#include "include/event/dssd_event.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run end_run\n"
		"  run               Set run number.\n"
		"  end_run           Set the last run to chain, included.\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -t tag            Set trigger tag.\n";
}

/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] trigger_tag trigger tag get from arguments
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	std::string &trigger_tag
) {
	// initialize
	help = false;
	trigger_tag.clear();
	// start index of positional arugments
	int result = 0;
	for (result = 1; result < argc; ++result) {
		// assumed that all options have read
		if (argv[result][0] != '-') break;
		// short option contains only one letter
		if (argv[result][2] != 0) return -result;
		if (argv[result][1] == 'h') {
			help = true;
			return result;
		} else if (argv[result][1] == 't') {
			// option of trigger tag
			// get tag in next argument
			++result;
			// miss arguemnt behind option
			if (result == argc) return -argc;
			trigger_tag = argv[result];
		} else {
			return -result;
		}
	}
	return result;
}

int main(int argc, char **argv) {
	if (argc < 3) {
		PrintUsage(argv[0]);
		return -1;
	}
	// help flag
	bool help = false;
	// trigger tag
	std::string tag;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, tag);
	// need help
	if (help) {
		PrintUsage(argv[0]);
		return 0;
	}
	// invalid arguments
	if (pos_start < 0) {
		if (-pos_start < argc) {
			std::cerr << "Error: Invaild option " << argv[-pos_start] << ".\n";
		} else {
			std::cerr << "Error: Option need parameter.\n";
		}
		PrintUsage(argv[0]);
		return -1;
	}
	// check number of positional arguments
	if (pos_start+1 >= argc) {
		// positional arguments less than 3
		std::cerr << "Error: Miss detector argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}
	// run number
	unsigned int run = atoi(argv[pos_start]);
	// run length
	unsigned int end_run = atoi(argv[pos_start+1]);

	// input time reference chain
	TChain ref_chain("ref", "chain");
	for (unsigned int i = run; i <= end_run; ++i) {
		if (i == 628) continue;
		ref_chain.AddFile(TString::Format(
			"%s%sreftime-%s%04u.root/tree",
			kGenerateDataPath,
			kFundamentalDir,
			tag.empty() ? "" : (tag+"-").c_str(),
			i
		));
	}
	// input t0d1 chain
	TChain t0d1_chain("t0d1", "t0d1");
	for (unsigned int i = run; i <= end_run; ++i) {
		if (i == 628) continue;
		t0d1_chain.AddFile(TString::Format(
			"%s%st0d1-fundamental-%s%04u.root/tree",
			kGenerateDataPath,
			kFundamentalDir,
			tag.empty() ? "" : (tag+"-").c_str(),
			i
		));
	}
	// input t0d2 chain
	TChain t0d2_chain("t0d2", "t0d2");
	for (unsigned int i = run; i <= end_run; ++i) {
		if (i == 628) continue;
		t0d2_chain.AddFile(TString::Format(
			"%s%st0d2-fundamental-%s%04u.root/tree",
			kGenerateDataPath,
			kFundamentalDir,
			tag.empty() ? "" : (tag+"-").c_str(),
			i
		));
	}
	// input t0d3 chain
	TChain t0d3_chain("t0d3", "t0d3");
	for (unsigned int i = run; i <= end_run; ++i) {
		if (i == 628) continue;
		t0d3_chain.AddFile(TString::Format(
			"%s%st0d3-fundamental-%s%04u.root/tree",
			kGenerateDataPath,
			kFundamentalDir,
			tag.empty() ? "" : (tag+"-").c_str(),
			i
		));
	}
	// add friend
	ref_chain.AddFriend(&t0d1_chain);
	ref_chain.AddFriend(&t0d2_chain);
	ref_chain.AddFriend(&t0d3_chain);
	// input reference time
	double ref_time;
	// input t0d1 fundamental event
	DssdFundamentalEvent t0d1;
	// input t0d2 fundamental event
	DssdFundamentalEvent t0d2;
	// input t0d3 fundamental event
	DssdFundamentalEvent t0d3;
	// setup input branches
	ref_chain.SetBranchAddress("time", &ref_time);
	t0d1.SetupInput(&ref_chain, "t0d1.");
	t0d2.SetupInput(&ref_chain, "t0d2.");
	t0d3.SetupInput(&ref_chain, "t0d3.");

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%stime-energy-%s%04u-%04u.root",
		kGenerateDataPath,
		kShowDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run,
		end_run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// time-energy histogram of t0d1 front strips [0, 32), [48, 64)
	TH2F t0d1f1(
		"hd1f1", "t0d1 front strips [0, 32), [48, 64)",
		1000, 0, 60000, 1000, -200, 800
	);
	// time-energy histogram of t0d1 front strips [32, 48)
	TH2F t0d1f2(
		"hd1f2", "t0d1 front strips [32, 48)",
		1000, 0, 60000, 1000, -200, 800
	);
	// time-energy histogram of t0d1 back strips
	TH2F t0d1b(
		"hd1b", "t0d1 back strips",
		1000, 0, 60000, 1000, -200, 800
	);
	// time-energy histogram of t0d2 front strips
	TH2F t0d2f(
		"hd2f", "t0d2 front strips",
		1000, 0, 60000, 1000, -200, 800
	);
	// time-energy histogram of t0d2 back strips
	TH2F t0d2b(
		"hd2b", "t0d2 back strips",
		1000, 0, 60000, 1000, -200, 800
	);
	// time-energy histogram of t0d3 front strips
	TH2F t0d3f(
		"hd3f", "t0d3 front strips",
		1000, 0, 60000, 1000, -200, 800
	);
	// time-energy histogram of t0d3 back strips
	TH2F t0d3b(
		"hd3b", "t0d3 back strips",
		1000, 0, 60000, 1000, -200, 800
	);

	// total number of entries
	long long entries = ref_chain.GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling time-energy histogram   0%%");
	fflush(stdout);
	// loop to fill histogram
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ref_chain.GetEntry(entry);
		if (ref_time < -9e4) continue;
		// fill t0d1
		for (unsigned short i = 0; i < t0d1.front_hit; ++i) {
			if (t0d1.front_strip[i] < 32 || t0d1.front_strip[i] >= 48) {
				t0d1f1.Fill(t0d1.front_energy[i], t0d1.front_time[i]-ref_time);
			} else {
				t0d1f2.Fill(t0d1.front_energy[i], t0d1.front_time[i]-ref_time);
			}
		}
		for (unsigned short i = 0; i < t0d1.back_hit; ++i) {
			t0d1b.Fill(t0d1.back_energy[i], t0d1.back_time[i]-ref_time);
		}
		// fill t0d2
		for (unsigned short i = 0; i < t0d2.front_hit; ++i) {
			t0d2f.Fill(t0d2.front_energy[i], t0d2.front_time[i]-ref_time);
		}
		for (unsigned short i = 0; i < t0d2.back_hit; ++i) {
			t0d2b.Fill(t0d2.back_energy[i], t0d2.back_time[i]-ref_time);
		}
		// fill t0d3
		for (unsigned short i = 0; i < t0d3.front_hit; ++i) {
			t0d3f.Fill(t0d3.front_energy[i], t0d3.front_time[i]-ref_time);
		}
		for (unsigned short i = 0; i < t0d3.back_hit; ++i) {
			t0d3b.Fill(t0d3.back_energy[i], t0d3.back_time[i]-ref_time);
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save time-energy histograms
	t0d1f1.Write();
	t0d1f2.Write();
	t0d1b.Write();
	t0d2f.Write();
	t0d2b.Write();
	t0d3f.Write();
	t0d3b.Write();
	// close files
	opf.Close();
	return 0;
}