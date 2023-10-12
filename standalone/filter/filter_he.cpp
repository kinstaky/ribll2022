#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TString.h>

#include "include/event/dssd_event.h"
#include "include/event/t0_event.h"
#include "include/event/particle_type_event.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run detector\n"
		"  run               Set run number.\n"
		"  detector          Set detector name.\n"
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
		// positional arguments less than 1
		std::cerr << "Error: Miss detector argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}
	// run number
	unsigned int run = atoi(argv[pos_start]);
	// detector name
	std::string detector(argv[pos_start+1]);
	if (detector != "t0d1" && detector != "t0d2" && detector != "t0d3") {
		return -1;
	}
	unsigned short layer = argv[pos_start+1][3] - '0';

	// input telescope file name
	TString t0_file_name;
	t0_file_name.Form(
		"%s%st0-telescope-%s%04u.root",
		kGenerateDataPath,
		kTelescopeDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// input file
	TFile t0_file(t0_file_name, "read");
	// input t0 tree
	TTree *ipt = (TTree*)t0_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< t0_file_name << " failed.\n";
		return -1;
	}
	// add particle type friends
	ipt->AddFriend(
		"type=tree",
		TString::Format(
			"%s%st0-particle-type-%s%04u.root",
			kGenerateDataPath,
			kParticleIdentifyDir,
			tag.empty() ? "" : (tag+"-").c_str(),
			run
		)
	);
	// add dssd normalize result friends
	ipt->AddFriend(
		"result=tree",
		TString::Format(
			"%s%s%s-result-%s%04u.root",
			kGenerateDataPath,
			kNormalizeDir,
			detector.c_str(),
			tag.empty() ? "" : (tag+"-").c_str(),
			run
		)
	);
	// add dssd fundamental friends
	ipt->AddFriend(
		"fundamental=tree",
		TString::Format(
			"%s%s%s-fundamental-%s%04u.root",
			kGenerateDataPath,
			kFundamentalDir,
			detector.c_str(),
			tag.empty() ? "" : (tag+"-").c_str(),
			run
		)
	);
	// input T0 event
	T0Event t0_event;
	// input particle type event
	ParticleTypeEvent type_event;
	// DSSD normalize result event
	DssdFundamentalEvent result_event;
	// DSSD fundamental event
	DssdFundamentalEvent fundamental_event;
	// setup input branches
	t0_event.SetupInput(ipt);
	type_event.SetupInput(ipt, "type.");
	result_event.SetupInput(ipt, "result.");
	fundamental_event.SetupInput(ipt, "fundamental.");

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-fundamental-%she-%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		detector.c_str(),
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "filter He events");
	// output event
	DssdFundamentalEvent output_event;
	// setup output branches
	output_event.SetupOutput(&opt);

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filtering events   0%%");
	fflush(stdout);
	// loop to fill pid histogram
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);
		for (unsigned short i = 0; i < t0_event.num; ++i) {
			if (type_event.charge[i] != 2) continue;
			if (type_event.mass[i] != 4) continue;
			if (t0_event.layer[i] != 2) continue;
			if (((t0_event.status[i] >> (4 * (layer-1))) & 0xf) != 0) continue;
			unsigned int front_flag = 0;
			unsigned int back_flag = 0;
			if (layer == 1) {
				front_flag = t0_event.dssd_flag[i][0] & 0xff;
				back_flag = (t0_event.dssd_flag[i][0] >> 8) & 0xff;
			} else if (layer == 2) {
				front_flag = t0_event.dssd_flag[i][1] & 0xff;
				back_flag = (t0_event.dssd_flag[i][1] >> 8) & 0xff;
			} else if (layer == 3) {
				front_flag = t0_event.dssd_flag[i][2] & 0xff;
				back_flag = (t0_event.dssd_flag[i][2] >> 8) & 0xff;
			} else {
				std::cerr << "Error: Invalid Dssd layer " << layer << ".\n";
				return -1;
			}
			int result_front_index = -1;
			for (unsigned int i = 0; i < 8; ++i) {
				if ((front_flag >> i) == 1) {
					result_front_index = i;
					break;
				}
			}
			int result_back_index = -1;
			for (unsigned int i = 0; i < 8; ++i) {
				if ((back_flag >> i) == 1) {
					result_back_index = i;
					break;
				}
			}
			if (result_front_index < 0) continue;
			if (result_front_index > result_event.front_hit) continue;
			if (result_back_index < 0) continue;
			if (result_back_index > result_event.back_hit) continue;

			// fill event
			output_event.front_hit = 1;
			output_event.back_hit = 1;
			output_event.cfd_flag = 0;
			unsigned short front_index =
				result_event.front_fundamental_index[result_front_index];
			unsigned short back_index =
				result_event.back_fundamental_index[result_back_index];
			output_event.cfd_flag |=
				(fundamental_event.cfd_flag >> front_index) & 0x1;
			output_event.cfd_flag |=
				(fundamental_event.cfd_flag >> (back_index)) & 0x100;
			output_event.front_strip[0]
				= fundamental_event.front_strip[front_index];
			output_event.back_strip[0]
				= fundamental_event.back_strip[back_index];
			output_event.front_time[0]
				= fundamental_event.front_time[front_index];
			output_event.back_time[0]
				= fundamental_event.back_time[back_index];
			output_event.front_energy[0]
				= fundamental_event.front_energy[front_index];
			output_event.back_energy[0]
				= fundamental_event.back_energy[back_index];
			output_event.front_decode_entry[0]
				= fundamental_event.front_decode_entry[front_index];
			output_event.back_decode_entry[0]
				= fundamental_event.back_decode_entry[back_index];
// if (entry < 100000) {
// 	std::cout << entry << " " << output_event.front_energy[0] << " " << output_event.back_energy[0] << "\n";
// }
			opt.Fill();
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save tree
	opt.Write();
	// close files
	opf.Close();
	return 0;
}