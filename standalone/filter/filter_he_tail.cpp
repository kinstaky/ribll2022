#include <iostream>
#include <vector>
#include <fstream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>
#include <TCutG.h>

#include "include/event/dssd_event.h"
#include "include/event/filter_event.h"
#include "include/event/particle_type_event.h"
#include "include/event/t0_event.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run\n"
		"  run               Set run number.\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -t tag            Set trigger tag.\n"
		"  -i num            Set iteration mode.\n";
}

/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] trigger_tag trigger tag get from arguments
/// @param[out] iteartion iteration mode
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	std::string &trigger_tag,
	int &iteration
) {
	// initialize
	help = false;
	trigger_tag.clear();
	iteration = 0;
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
		} else if (argv[result][1] == 'i') {
			// option of iteration flag
			// get number in next argument
			++result;
			// miss arguemnt behind option
			if (result == argc) return -argc;
			iteration = atoi(argv[result]);
		} else {
			return -result;
		}
	}
	return result;
}


/// @brief search for index from binary flag
/// @param[in] side 0-front, 1-back
/// @param[in] flag binary flag to search
/// @returns valid index of front or back side if only one valid bit is found,
/// 	otherwise, returns -1
int SearchIndex(int side, unsigned short flag) {
	int counts = 0;
	int index = -1;
	int base = side == 0 ? 0x1 : 0x100;
	for (int i = 0; i < 8; ++i) {
		if ((flag & (base << i)) != 0) {
			++counts;
			index = i;
		}
	}
	if (counts != 1) return -1;
	return index;
}


std::unique_ptr<TCutG> ReadCut(const char *name) {

	// cut file name
	TString cut_file_name;
	cut_file_name.Form(
		"%s%scut/%s.txt",
		kGenerateDataPath,
		kParticleIdentifyDir,
		name
	);
	// open cut file to read points
	std::ifstream fin(cut_file_name.Data());
	if (!fin.good()) {
		std::cerr << "Error: Open file "
			<< cut_file_name << " failed.\n";
		return nullptr;
	}
	// result
	std::unique_ptr<TCutG> result = std::make_unique<TCutG>();
	// point index
	int point;
	// point positions
	double x, y;
	// loop to read points
	while (fin.good()) {
		fin >> point >> x >> y;
		result->SetPoint(point, x, y);
	}
	// close file
	fin.close();
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
	// iteration flag
	int iteration = 0;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, tag, iteration);

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
		// positional arguments less than 3
		std::cerr << "Error: Miss detector argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}

	int run = atoi(argv[pos_start]);

	// input T0 telescope file name
	TString t0_telescope_file_name = TString::Format(
		"%s%st0-telescope-%s%04d.root",
		kGenerateDataPath,
		kTelescopeDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// input T0 particle type file
	TFile t0_telescope_file(t0_telescope_file_name, "read");
	// input T0 particle tree
	TTree *ipt = (TTree*)t0_telescope_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< t0_telescope_file_name << " failed.\n";
		return -1;
	}

	// loop to add T0Dx friends
	for (int i = 1; i <= 3; ++i) {
		// input T0Dx normalize result file name
		TString norm_result_file_name = TString::Format(
			"%s%st0d%d-result-%s%04d.root",
			kGenerateDataPath,
			kNormalizeDir,
			i,
			tag.empty() ? "" : (tag+"-").c_str(),
			run
		);
		// add T0Dx normalize result friend
		ipt->AddFriend(
			TString::Format("d%dr=tree", i),
			norm_result_file_name
		);

		// // input T0Dx fundamental file name
		// TString fundamental_file_name = TString::Format(
		// 	"%s%st0d%d-fundamental-%s%04d.root",
		// 	kGenerateDataPath,
		// 	kFundamentalDir,
		// 	i,
		// 	tag.empty() ? "" : (tag+"-").c_str(),
		// 	run
		// );
		// // add T0Dx fundamental friend
		// ipt->AddFriend(
		// 	TString::Format("d%df=tree", i),
		// 	fundamental_file_name
		// );
	}

	// input t0 event
	T0Event t0_event;
	// input T0Dx normalize result events
	DssdFundamentalEvent norm_events[3];
	// // input T0Dx fundamental events
	// DssdFundamentalEvent fundamental_events[3];

	// setup input branches
	t0_event.SetupInput(ipt);
	for (int i = 0; i < 3; ++i) {
		norm_events[i].SetupInput(
			ipt, TString::Format("d%dr.", i+1).Data()
		);
		// fundamental_events[i].SetupInput(
		// 	ipt, TString::Format("d%df.", i+1).Data()
		// );
	}

	TString output_file_name = TString::Format(
		"%s%st0-he-tail-%s%04d.root",
		kGenerateDataPath,
		kFilterDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// output file
	TFile output_file(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "4He tail information");
	// output data
	unsigned short num;
	unsigned short flag[4];
	double energy[4][3];
	double x[4][3];
	double y[4][3];
	unsigned short merge_flag[4][3];
	unsigned short front_hit[3];
	unsigned short back_hit[3];
	unsigned short front_strip[3][8];
	unsigned short back_strip[3][8];
	double front_energy[3][8];
	double back_energy[3][8];
	double front_time[3][8];
	double back_time[3][8];
	// setup output branches
	opt.Branch("num", &num, "num/s");
	opt.Branch("flag", flag, "flag[num]/s");
	opt.Branch("energy", energy, "e[num][3]/D");
	opt.Branch("x", x, "x[num][3]/D");
	opt.Branch("y", y, "y[num][3]/D");
	opt.Branch("merge_flag", merge_flag, "mflag[num][3]/s");
	for (int i = 0; i < 3; ++i) {
		opt.Branch(
			TString::Format("d%dfhit", i+1),
			front_hit+i,
			TString::Format("d%dfhit/s", i+1)
		);
		opt.Branch(
			TString::Format("d%dbhit", i+1),
			back_hit+i,
			TString::Format("d%dbhit/s", i+1)
		);
		opt.Branch(
			TString::Format("d%dfs", i+1),
			front_strip[i],
			TString::Format("d%dfs[d%dfhit]/s", i+1, i+1)
		);
		opt.Branch(
			TString::Format("d%dbs", i+1),
			back_strip[i],
			TString::Format("d%dbs[d%dbhit]/s", i+1, i+1)
		);
		opt.Branch(
			TString::Format("d%dfe", i+1),
			front_energy[i],
			TString::Format("d%dfe[d%dfhit]/D", i+1, i+1)
		);
		opt.Branch(
			TString::Format("d%dbe", i+1),
			back_energy[i],
			TString::Format("d%dbe[d%dbhit]/D", i+1, i+1)
		);
		opt.Branch(
			TString::Format("d%dft", i+1),
			front_time[i],
			TString::Format("d%dft[d%dfhit]/D", i+1, i+1)
		);
		opt.Branch(
			TString::Format("d%dbt", i+1),
			back_time[i],
			TString::Format("d%dbt[d%dbhit]/D", i+1, i+1)
		);
	}

	// read He tail cut
	std::unique_ptr<TCutG> tail_cut = ReadCut("t0-d1d2-tail-p1i1-He");

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filtering identifying events   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get events
		ipt->GetEntry(entry);

		if (t0_event.num != 1 && t0_event.num != 2) continue;
		num = t0_event.num;
		// check is invalid tail
		bool is_tail = false;
		for (unsigned short i = 0; i < num; ++i) {
			if (
				t0_event.flag[i] == 0x3
				&& tail_cut->IsInside(
					t0_event.energy[i][1], t0_event.energy[i][0]
				)
			) {
				is_tail = true;
				break;
			}
		}
		if (!is_tail) continue;
		// rec
		for (unsigned short i = 0; i < num; ++i) {
			flag[i] = t0_event.flag[i];
			for (int j = 0; j < 3; ++j) {
				energy[i][j] = t0_event.energy[i][j];
				x[i][j] = t0_event.x[i][j];
				y[i][j] = t0_event.y[i][j];
			}
			merge_flag[i][0] = t0_event.d1_flag[i];
			merge_flag[i][1] = t0_event.d2_flag[i];
			merge_flag[i][2] = t0_event.d3_flag[i];
		}
		for (int i = 0; i < 3; ++i) {
			front_hit[i] = norm_events[i].front_hit;
			for (unsigned short j = 0; j < front_hit[i]; ++j) {
				front_strip[i][j] = norm_events[i].front_strip[j];
				front_energy[i][j] = norm_events[i].front_energy[j];
				front_time[i][j] = norm_events[i].front_time[j];
			}
			back_hit[i] = norm_events[i].back_hit;
			for (unsigned short j = 0; j < back_hit[i]; ++j) {
				back_strip[i][j] = norm_events[i].back_strip[j];
				back_energy[i][j] = norm_events[i].back_energy[j];
				back_time[i][j] = norm_events[i].back_time[j];
			}
		}
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save tree
	opt.Write();
	// close files
	t0_telescope_file.Close();
	output_file.Close();

	return 0;
}