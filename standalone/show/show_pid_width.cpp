#include <iostream>
#include <iomanip>
#include <fstream>

#include <TChain.h>
#include <TCutG.h>
#include <TFile.h>
#include <TGraph.h>
#include <TString.h>
#include <TTree.h>

#include "include/event/t0_event.h"
#include "include/event/particle_event.h"

using namespace ribll;


/// @brief print usage of this program
/// @param[in] name program name
///
void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run [end_run]\n"
		"  run               Set run number.\n"
		"  end_run           Set the last run, inclusive.\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -t tag            Set trigger tag.\n"
		"Examples:\n"
		"  " << name << " 600        Show center of run 600.\n"
		"  " << name << " 600 700    Show center from run 600 to 700.\n";
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
	if (argc < 2) {
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
	if (pos_start < 0) {
		if (-pos_start < argc) {
			std::cerr << "Error: Invaild option " << argv[-pos_start] << ".\n";
		} else {
			std::cerr << "Error: Option need parameter.\n";
		}
		PrintUsage(argv[0]);
		return -1;
	}
	if (pos_start >= argc) {
		// positional arguments less than 1
		std::cerr << "Error: Miss run argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}

	// run number
	unsigned int run = atoi(argv[pos_start]);
	unsigned int end_run = run;
	if (pos_start + 1 < argc) {
		end_run = atoi(argv[pos_start+1]);
	}

	// T0 event chain
	TChain ipt("t0", "chain of T0 events");
	for (unsigned int i = run; i <= end_run; ++i) {
		if (i == 628) continue;
		ipt.AddFile(TString::Format(
			"%s%st0-telescope-%s%04u.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			tag.empty() ? "" : (tag+"-").c_str(),
			i
		));
	}
	// PPAC event chain
	TChain ppac_chain("ppac", "chain of PPAC events");
	for (unsigned int i = run; i <= end_run; ++i) {
		if (i == 628) continue;
		ppac_chain.AddFile(TString::Format(
			"%s%sxppac-particle-%s%04u.root/tree",
			kGenerateDataPath,
			kParticleDir,
			tag.empty() ? "" : (tag+"-").c_str(),
			i
		));
	}
	// add friend
	ipt.AddFriend(&ppac_chain);
	// input telescope event
	T0Event t0_event;
	// input ppac event
	ParticleEvent ppac_event;
	// setup input branches
	t0_event.SetupInput(&ipt);
	ppac_event.SetupInput(&ipt, "ppac.");

	// output calibration root file name
	TString output_file_name;
	output_file_name.Form(
		"%s%st0-pid-width-%s%04u-%04u.root",
		kGenerateDataPath,
		kShowDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run,
		end_run
	);
	// output calibration file
	TFile opf(output_file_name, "recreate");
	// output E-VS-dE graph
	TGraph e_vs_de;

	// T0D1D2 4He cut
	std::ifstream fin(TString::Format(
		"%s%scut/t0-d1d2-4He.txt",
		kGenerateDataPath,
		kParticleIdentifyDir
	).Data());
	if (!fin.good()) {
		std::cerr << "Error: Read cut failed.\n";
		return -1;
	}
	// cut
	TCutG cut;
	// point index
	int points;
	// point positions
	double cutx, cuty;
	// read points
	while (fin.good()) {
		fin >> points >> cutx >> cuty;
		cut.SetPoint(points, cutx, cuty);
	}
	// close file
	fin.close();


 	// total number of entries
	long long entries = ipt.GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling events to graph   0%%");
	fflush(stdout);
	// loop to fill events to graph
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt.GetEntry(entry);
		for (unsigned short i = 0; i < t0_event.num; ++i) {
			if (
				cut.IsInside(t0_event.energy[i][1], t0_event.energy[i][0])
				&& t0_event.status[i] == 1111
				&& fabs(t0_event.x[i][0] - t0_event.x[i][1]) < 2
				&& fabs(t0_event.y[i][0] - t0_event.y[i][1]) < 2
				&& ppac_event.num == 4
				&& fabs(ppac_event.x[3] - t0_event.x[i][0]) < 2
				&& fabs(ppac_event.y[3] - t0_event.y[i][0]) < 2
			) {
				e_vs_de.AddPoint(t0_event.energy[i][0], t0_event.energy[i][1]);
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save graph
	opf.cd();
	e_vs_de.Write("gede");
	// close files
	opf.Close();

	return 0;
}
