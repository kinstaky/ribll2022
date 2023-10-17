#include <iostream>
#include <vector>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2F.h>

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

	// input T0 particle type file name
	TString pid_file_name = TString::Format(
		"%s%st0-particle-type-%s%04d.root",
		kGenerateDataPath,
		kParticleIdentifyDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// input T0 particle type file
	TFile pid_file(pid_file_name, "read");
	// input T0 particle tree
	TTree *ipt = (TTree*)pid_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< pid_file_name << " failed.\n";
		return -1;
	}

	// input T0 telescope file name
	TString t0_telescope_file_name = TString::Format(
		"%s%st0-telescope-%s%04d.root",
		kGenerateDataPath,
		kTelescopeDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// add T0 telescope friend
	ipt->AddFriend("t0=tree", t0_telescope_file_name);

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

		// input T0Dx fundamental file name
		TString fundamental_file_name = TString::Format(
			"%s%st0d%d-fundamental-%s%04d.root",
			kGenerateDataPath,
			kFundamentalDir,
			i,
			tag.empty() ? "" : (tag+"-").c_str(),
			run
		);
		// add T0Dx fundamental friend
		ipt->AddFriend(
			TString::Format("d%df=tree", i),
			fundamental_file_name
		);
	}

	// input particle type event
	ParticleTypeEvent pid_event;
	// input t0 event
	T0Event t0_event;
	// input T0Dx normalize result events
	DssdFundamentalEvent norm_events[3];
	// input T0Dx fundamental events
	DssdFundamentalEvent fundamental_events[3];

	// setup input branches
	pid_event.SetupInput(ipt);
	t0_event.SetupInput(ipt, "t0.");
	for (int i = 0; i < 3; ++i) {
		norm_events[i].SetupInput(
			ipt, TString::Format("d%dr.", i+1).Data()
		);
		fundamental_events[i].SetupInput(
			ipt, TString::Format("d%df.", i+1).Data()
		);
	}

	// output filter files
	TFile* output_files[3];
	// histogram of strips
	TH2F* hist_strip[3];
	// histogram of energy before normalization
	TH2F* hist_energy[3];
	// output filter trees
	TTree* opts[3];
	// output filter events
	FilterEvent filter_events[3];
	for (int i = 0; i < 3; ++i) {
		TString output_file_name = TString::Format(
			"%s%st0d%d-normalize-filter-1-%s%04d.root",
			kGenerateDataPath,
			kFilterDir,
			i+1,
			tag.empty() ? "" : (tag+"-").c_str(),
			run
		);
		// output file
		output_files[i] = new TFile(output_file_name, "recreate");
		int strips = i == 0 ? 64 : 32;
		// histograms
		hist_strip[i] = new TH2F(
			"hfbs", "strips", strips, 0, strips, strips, 0, strips
		);
		hist_energy[i] = new TH2F(
			"hfbe", "front back energy before normalization",
			1000, 0, 60000,
			1000, 0, 60000
		);
		// output tree
		opts[i] = new TTree(
			"tree", TString::Format("T0D%d normalize filter 1", i+1)
		);
		// setup output branches
		filter_events[i].SetupOutput(opts[i]);
	}

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

		// initialize filter events
		for (int i = 0; i < 3; ++i) filter_events[i].num = 0;

		for (unsigned short i = 0; i < pid_event.num; ++i) {
			if (pid_event.layer[i] < 1) continue;
			if (pid_event.mass[i] == 0 || pid_event.charge[i] == 0) continue;

			for (int j = 0; j < (pid_event.layer[i] >= 2 ? 3 : 2); ++j) {
				// for convenience
				int &num = filter_events[j].num;

				// index in particle type event
				filter_events[j].pid_index[num] = i;
				// merge flag
				filter_events[j].merge_flag[num] = t0_event.dssd_flag[i][j];

				// get normalize result event index from merge flag
				int front_norm_index =
					SearchIndex(0, t0_event.dssd_flag[i][j]);
				int back_norm_index =
					SearchIndex(1, t0_event.dssd_flag[i][j]);
				// check index
				if (front_norm_index < 0 || back_norm_index < 0) continue;
				// valid index, fill to filter event
				filter_events[j].norm_front_index[num] = front_norm_index;
				filter_events[j].norm_back_index[num] = back_norm_index;

				// get fundamental event index from normalize result event
				unsigned short front_fundamental_index =
					norm_events[j].front_fundamental_index[front_norm_index];
				unsigned short back_fundamental_index =
					norm_events[j].back_fundamental_index[back_norm_index];
				// fill to filter event
				filter_events[j].front_index[num] = front_fundamental_index;
				filter_events[j].back_index[num] = back_fundamental_index;

				// get strips from fundamental events
				unsigned short front_strip =
					fundamental_events[j].
						front_strip[front_fundamental_index];
				unsigned short back_strip =
					fundamental_events[j].
						back_strip[back_fundamental_index];
				// get energy from fundamental events
				double front_energy =
					fundamental_events[j].
						front_energy[front_fundamental_index];
				double back_energy =
					fundamental_events[j].
						back_energy[back_fundamental_index];
				// fill to filter events
				filter_events[j].front_strip[num] = front_strip;
				filter_events[j].back_strip[num] = back_strip;
				filter_events[j].front_energy[num] = front_energy;
				filter_events[j].back_energy[num] = back_energy;

				// fill to histogram
				hist_strip[j]->Fill(front_strip, back_strip);
				hist_energy[j]->Fill(front_energy, back_energy);

				// update num
				++filter_events[j].num;
			}
		}

		// fill filter event
		for (int i = 0; i < 3; ++i) {
			output_files[i]->cd();
			opts[i]->Fill();
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save histograms and trees
	for (int i = 0; i < 3; ++i) {
		output_files[i]->cd();
		hist_strip[i]->Write();
		hist_energy[i]->Write();
		opts[i]->Write();
	}

	// close files
	pid_file.Close();
	for (auto &file : output_files) {
		file->Close();
	}

	return 0;
}