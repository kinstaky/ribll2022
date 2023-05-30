#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TString.h>

#include "include/event/dssd_event.h"
#include "include/detectors.h"

using namespace ribll;

void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run end_run detector\n"
		"  run               Set run number.\n"
		"  end_run           Set the last run to chain, included.\n"
		"  detector          Set the detector name.\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -t tag            Set trigger tag.\n"
		"  -n                Use normalized data.\n";
}

/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] trigger_tag trigger tag get from arguments
/// @param[out] normalize use normalized data
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	std::string &trigger_tag,
	bool &normalize
) {
	// initialize
	help = false;
	trigger_tag.clear();
	normalize = false;
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
		} else if (argv[result][1] == 'n') {
			normalize = true;
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
	// normalize
	bool normalize = false;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, tag, normalize);
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
	if (pos_start+2 >= argc) {
		// positional arguments less than 3
		std::cerr << "Error: Miss detector argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}
	// run number
	unsigned int run = atoi(argv[pos_start]);
	// run length
	unsigned int end_run = atoi(argv[pos_start+1]);
	// detector name
	std::string detector_name(argv[pos_start+2]);

	std::shared_ptr<Dssd> dssd = CreateDssd(detector_name, run, tag);

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
	// input detector chain
	TChain detector_chain(detector_name.c_str(), "detector");
	for (unsigned int i = run; i <= end_run; ++i) {
		if (i == 628) continue;
		if (!normalize) {
			detector_chain.AddFile(TString::Format(
				"%s%s%s-fundamental-%s%04u.root/tree",
				kGenerateDataPath,
				kFundamentalDir,
				detector_name.c_str(),
				tag.empty() ? "" : (tag+"-").c_str(),
				i
			));
		} else {
			detector_chain.AddFile(TString::Format(
				"%s%s%s-result-%s%04u-0.root/tree",
				kGenerateDataPath,
				kNormalizeDir,
				detector_name.c_str(),
				tag.empty() ? "" : (tag+"-").c_str(),
				i
			));
		}
	}
	// add friend
	ref_chain.AddFriend(&detector_chain);
	// input reference time
	double ref_time;
	// input dssd fundamental event
	DssdFundamentalEvent event;
	// setup input branches
	ref_chain.SetBranchAddress("time", &ref_time);
	event.SetupInput(&ref_chain, (detector_name+".").c_str());

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-time-energy-hit-%s%s%04u-%04u.root",
		kGenerateDataPath,
		kShowDir,
		detector_name.c_str(),
		tag.empty() ? "" : (tag+"-").c_str(),
		normalize ? "norm-" : "",
		run,
		end_run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// time-energy 2D histogram of single strip and single hit
	std::vector<TH2F> time_energy_single;
	for (size_t side = 0;  side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_energy_single.emplace_back(
				TString::Format("%cte%ld", "fb"[side], i),
				"time energy",
				1000, 0, 60000,
				2000, normalize ? -500 : -200,  normalize ? 500 : 800
			);
		}
	}
	// time-energy 2D histogram of single strip and single hit with valid CFD
	std::vector<TH2F> time_energy_single_cfd;
	for (size_t side = 0;  side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_energy_single_cfd.emplace_back(
				TString::Format("%ctec%ld", "fb"[side], i),
				"time energy with CFD",
				1000, 0, 60000,
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time-energy 2D histogram of single strip and single hit with LE
	std::vector<TH2F> time_energy_single_le;
	for (size_t side = 0;  side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_energy_single_le.emplace_back(
				TString::Format("%ctel%ld", "fb"[side], i),
				"time energy with LE",
				1000, 0, 60000,
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time 1D histogram of single strip and single hit
	std::vector<TH1F> time_single;
	for (size_t side = 0; side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_single.emplace_back(
				TString::Format("%ct%ld", "fb"[side], i),
				"time",
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time 1D histogram of single strip and single hit with CFD
	std::vector<TH1F> time_single_cfd;
	for (size_t side = 0; side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_single_cfd.emplace_back(
				TString::Format("%ctc%ld", "fb"[side], i),
				"time with CFD",
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time 1D histogram of single strip and single hit with LE
	std::vector<TH1F> time_single_le;
	for (size_t side = 0; side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_single_le.emplace_back(
				TString::Format("%ctl%ld", "fb"[side], i),
				"time with LE",
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}

	// time-energy 2D histogram of single strip and double hit small
	std::vector<TH2F> time_energy_small;
	for (size_t side = 0;  side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_energy_small.emplace_back(
				TString::Format("%ctes%ld", "fb"[side], i),
				"time energy",
				1000, 0, 60000,
				2000, normalize ? -500 : -200,  normalize ? 500 : 800
			);
		}
	}
	// time-energy 2D histogram of single strip and double hit small with CFD
	std::vector<TH2F> time_energy_small_cfd;
	for (size_t side = 0;  side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_energy_small_cfd.emplace_back(
				TString::Format("%ctesc%ld", "fb"[side], i),
				"time energy with CFD",
				1000, 0, 60000,
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time-energy 2D histogram of single strip and double hit small with LE
	std::vector<TH2F> time_energy_small_le;
	for (size_t side = 0;  side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_energy_small_le.emplace_back(
				TString::Format("%ctesl%ld", "fb"[side], i),
				"time energy with LE",
				1000, 0, 60000,
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time 1D histogram of single strip and double hit small
	std::vector<TH1F> time_small;
	for (size_t side = 0; side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_small.emplace_back(
				TString::Format("%cts%ld", "fb"[side], i),
				"time",
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time 1D histogram of single strip and double hit small with CFD
	std::vector<TH1F> time_small_cfd;
	for (size_t side = 0; side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_small_cfd.emplace_back(
				TString::Format("%ctsc%ld", "fb"[side], i),
				"time with CFD",
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time 1D histogram of single strip and double hit small with LE
	std::vector<TH1F> time_small_le;
	for (size_t side = 0; side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_small_le.emplace_back(
				TString::Format("%ctsl%ld", "fb"[side], i),
				"time with LE",
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}

	// time-energy 2D histogram of single strip and double hit big
	std::vector<TH2F> time_energy_big;
	for (size_t side = 0;  side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_energy_big.emplace_back(
				TString::Format("%cteb%ld", "fb"[side], i),
				"time energy",
				1000, 0, 60000,
				2000, normalize ? -500 : -200,  normalize ? 500 : 800
			);
		}
	}
	// time-energy 2D histogram of single strip and double hit big with CFD
	std::vector<TH2F> time_energy_big_cfd;
	for (size_t side = 0;  side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_energy_big_cfd.emplace_back(
				TString::Format("%ctebc%ld", "fb"[side], i),
				"time energy with CFD",
				1000, 0, 60000,
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time-energy 2D histogram of single strip and double hit big with LE
	std::vector<TH2F> time_energy_big_le;
	for (size_t side = 0;  side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_energy_big_le.emplace_back(
				TString::Format("%ctebl%ld", "fb"[side], i),
				"time energy with LE",
				1000, 0, 60000,
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time 1D histogram of single strip and double hit big
	std::vector<TH1F> time_big;
	for (size_t side = 0; side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_big.emplace_back(
				TString::Format("%ctb%ld", "fb"[side], i),
				"time",
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time 1D histogram of single strip and double hit big with CFD
	std::vector<TH1F> time_big_cfd;
	for (size_t side = 0; side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_big_cfd.emplace_back(
				TString::Format("%ctbc%ld", "fb"[side], i),
				"time with CFD",
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time 1D histogram of single strip and double hit big with LE
	std::vector<TH1F> time_big_le;
	for (size_t side = 0; side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_big_le.emplace_back(
				TString::Format("%ctbl%ld", "fb"[side], i),
				"time with LE",
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}

	// time-energy 2D histogram of single module and single hit
	std::vector<TH2F> time_energy_single_module;
	for (size_t i = 0; i < 8; ++i) {
		time_energy_single_module.emplace_back(
			TString::Format("mte%ld", i),
			"time-energy of single module",
			1000, 0, 60000,
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time-energy 2D histogram of single module and single hit with CFD
	std::vector<TH2F> time_energy_single_module_cfd;
	for (size_t i = 0; i < 8; ++i) {
		time_energy_single_module_cfd.emplace_back(
			TString::Format("mtec%ld", i),
			"time-energy of single module with CFD",
			1000, 0, 60000,
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time-energy 2D histogram of single module and single hit with LE
	std::vector<TH2F> time_energy_single_module_le;
	for (size_t i = 0; i < 8; ++i) {
		time_energy_single_module_le.emplace_back(
			TString::Format("mtel%ld", i),
			"time-energy of single module with LE",
			1000, 0, 60000,
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time 1D histogram of single module and single hit
	std::vector<TH1F> time_single_module;
	for (size_t i = 0; i < 8; ++i) {
		time_single_module.emplace_back(
			TString::Format("mt%ld", i),
			"time of single module",
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time 1D histogram of single module and single hit with CFD
	std::vector<TH1F> time_single_module_cfd;
	for (size_t i = 0; i < 8; ++i) {
		time_single_module_cfd.emplace_back(
			TString::Format("mtc%ld", i),
			"time of single module with CFD",
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time 1D histogram of single module and single hit with LE
	std::vector<TH1F> time_single_module_le;
	for (size_t i = 0; i < 8; ++i) {
		time_single_module_le.emplace_back(
			TString::Format("mtl%ld", i),
			"time of single module with LE",
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}

	// time-energy 2D histogram of single module and double hit small
	std::vector<TH2F> time_energy_small_module;
	for (size_t i = 0; i < 8; ++i) {
		time_energy_small_module.emplace_back(
			TString::Format("mtes%ld", i),
			"time-energy of single module",
			1000, 0, 60000,
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time-energy 2D histogram of single module and doublt hit small with CFD
	std::vector<TH2F> time_energy_small_module_cfd;
	for (size_t i = 0; i < 8; ++i) {
		time_energy_small_module_cfd.emplace_back(
			TString::Format("mtesc%ld", i),
			"time-energy of single module with CFD",
			1000, 0, 60000,
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time-energy 2D histogram of single module and double hit small with LE
	std::vector<TH2F> time_energy_small_module_le;
	for (size_t i = 0; i < 8; ++i) {
		time_energy_small_module_le.emplace_back(
			TString::Format("mtesl%ld", i),
			"time-energy of single module with LE",
			1000, 0, 60000,
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time 1D histogram of single module and double hit small
	std::vector<TH1F> time_small_module;
	for (size_t i = 0; i < 8; ++i) {
		time_small_module.emplace_back(
			TString::Format("mts%ld", i),
			"time of single module",
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time 1D histogram of single module and double hit small with CFD
	std::vector<TH1F> time_small_module_cfd;
	for (size_t i = 0; i < 8; ++i) {
		time_small_module_cfd.emplace_back(
			TString::Format("mtsc%ld", i),
			"time of single module with CFD",
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time 1D histogram of single module and double hit small with LE
	std::vector<TH1F> time_small_module_le;
	for (size_t i = 0; i < 8; ++i) {
		time_small_module_le.emplace_back(
			TString::Format("mtsl%ld", i),
			"time of single module with LE",
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}

	// time-energy 2D histogram of single module and double hit big
	std::vector<TH2F> time_energy_big_module;
	for (size_t i = 0; i < 8; ++i) {
		time_energy_big_module.emplace_back(
			TString::Format("mteb%ld", i),
			"time-energy of single module",
			1000, 0, 60000,
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time-energy 2D histogram of single module and double hit big with CFD
	std::vector<TH2F> time_energy_big_module_cfd;
	for (size_t i = 0; i < 8; ++i) {
		time_energy_big_module_cfd.emplace_back(
			TString::Format("mtebc%ld", i),
			"time-energy of single module with CFD",
			1000, 0, 60000,
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time-energy 2D histogram of single module and double hit big with LE
	std::vector<TH2F> time_energy_big_module_le;
	for (size_t i = 0; i < 8; ++i) {
		time_energy_big_module_le.emplace_back(
			TString::Format("mtebl%ld", i),
			"time-energy of single module with LE",
			1000, 0, 60000,
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time 1D histogram of single module and double hit big
	std::vector<TH1F> time_big_module;
	for (size_t i = 0; i < 8; ++i) {
		time_big_module.emplace_back(
			TString::Format("mtb%ld", i),
			"time of single module",
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time 1D histogram of single module and double hit big with CFD
	std::vector<TH1F> time_big_module_cfd;
	for (size_t i = 0; i < 8; ++i) {
		time_big_module_cfd.emplace_back(
			TString::Format("mtbc%ld", i),
			"time of single module with CFD",
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time 1D histogram of single module and double hit big with LE
	std::vector<TH1F> time_big_module_le;
	for (size_t i = 0; i < 8; ++i) {
		time_big_module_le.emplace_back(
			TString::Format("mtbl%ld", i),
			"time of single module with LE",
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}

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
		if (
			event.front_hit == 1
			&& event.back_hit == 1
			&& (
				!normalize
				||
				fabs(event.front_energy[0]-event.back_energy[0]) < 500
			)
		) {
			unsigned short fs = event.front_strip[0];
			unsigned short bs = event.back_strip[0] + dssd->Strip(0);
			double fe = event.front_energy[0];
			double be = event.back_energy[0];
			double ft = event.front_time[0] - ref_time;
			double bt = event.back_time[0] - ref_time;
			unsigned short fm = fs / 16;
			unsigned short bm = bs / 16;
			// fill single strip histograms
			time_energy_single[fs].Fill(fe, ft);
			time_energy_single[bs].Fill(be, bt);
			time_single[fs].Fill(ft);
			time_single[bs].Fill(bt);
			// fill single module histograms
			time_energy_single_module[fm].Fill(fe, ft);
			time_energy_single_module[bm].Fill(be, bt);
			time_single_module[fm].Fill(ft);
			time_single_module[bm].Fill(bt);
			// check front strip cfd flag
			if ((event.cfd_flag & 0x1) == 0) {
				// fill single strip histograms
				time_energy_single_cfd[fs].Fill(fe, ft);
				time_single_cfd[fs].Fill(ft);
				// fill single module histograms
				time_energy_single_module_cfd[fm].Fill(fe, ft);
				time_single_module_cfd[fm].Fill(ft);
			} else {
				// fill single strip histograms
				time_energy_single_le[fs].Fill(fe, ft);
				time_single_le[fs].Fill(ft);
				// fill single module histograms
				time_energy_single_module_le[fm].Fill(fe, ft);
				time_single_module_le[fm].Fill(ft);
			}
			// check back strip cfd flag
			if ((event.cfd_flag & 0x100) == 0) {
				// fill single strip histograms
				time_energy_single_cfd[bs].Fill(be, bt);
				time_single_cfd[bs].Fill(bt);
				// fill single module histograms
				time_energy_single_module_cfd[bm].Fill(be, bt);
				time_single_module_cfd[bm].Fill(bt);
			} else {
				// fill single strip histograms
				time_energy_single_le[bs].Fill(be, bt);
				time_single_le[bs].Fill(bt);
				// fill single module histograms
				time_energy_single_module_le[bm].Fill(be, bt);
				time_single_module_le[bm].Fill(bt);
			}

		} else if (
			event.front_hit == 2
			&& event.back_hit == 1
			&& abs(event.front_strip[0]-event.front_strip[1]) == 1
			&& (
				!normalize
				||
				fabs(
					event.front_energy[0] + event.front_energy[1]
					- event.back_energy[0]
				) < 500
			)
		) {
			// small strip index
			size_t small = 0;
			// big strip index
			size_t big = 1;
			// swap if necessary
			if (event.front_strip[small] > event.front_strip[big]) {
				small = 1;
				big = 0;
			}
			// small strip energy
			double se = event.front_energy[small];
			// big strip energy
			double be = event.front_energy[big];
			// small strip time
			double st = event.front_time[small] - ref_time;
			// big strip time
			double bt = event.front_time[big] - ref_time;
			// small strip
			unsigned short ss = event.front_strip[small];
			// big strip
			unsigned short bs = event.front_strip[big];
			// small strip module
			size_t smod = ss == dssd->Strip(0)-1 ?
				ss / 16
				: (ss+1) / 16;
			// big strip module
			size_t bmod = bs == 0 ? 0 : (bs-1) / 16;

			// fill single strip histograms
			time_energy_small[ss].Fill(se, st);
			time_small[ss].Fill(st);
			time_energy_big[bs].Fill(be, bt);
			time_big[bs].Fill(bt);
			// fill single module histograms
			time_energy_small_module[smod].Fill(se, st);
			time_small_module[smod].Fill(st);
			time_energy_big_module[bmod].Fill(be, bt);
			time_big_module[bmod].Fill(bt);
			// check cfd flag
			if ((event.cfd_flag & (1 << small)) == 0) {
				// fill single strip histograms
				time_energy_small_cfd[ss].Fill(se, st);
				time_small_cfd[ss].Fill(st);
				// fill single module histograms
				time_energy_small_module_cfd[smod].Fill(se, st);
				time_small_module_cfd[smod].Fill(st);
			} else {
				// fill single strip histograms
				time_energy_small_le[ss].Fill(se, st);
				time_small_le[ss].Fill(st);
				// fill single module histograms
				time_energy_small_module_le[smod].Fill(se, st);
				time_small_module_le[smod].Fill(st);
			}
			if ((event.cfd_flag & (1 << big)) == 0) {
				// fill single strip histograms
				time_energy_big_cfd[bs].Fill(be, bt);
				time_big_cfd[bs].Fill(bt);
				// fill single module histograms
				time_energy_big_module_cfd[bmod].Fill(be, bt);
				time_big_module_cfd[bmod].Fill(bt);
			} else {
				// fill single strip histograms
				time_energy_big_le[bs].Fill(be, bt);
				time_big_le[bs].Fill(bt);
				// fill single module histograms
				time_energy_big_module_le[bmod].Fill(be, bt);
				time_big_module_le[bmod].Fill(bt);
			}

		} else if (
			event.front_hit == 1
			&& event.back_hit == 2
			&& abs(event.back_strip[0]-event.back_strip[1]) == 1
			&& (
				!normalize
				||
				fabs(
					event.front_energy[0]
					- event.back_energy[0] - event.back_energy[1]
				) < 500
			)
		) {
			// small strip index
			size_t small = 0;
			// big strip index
			size_t big = 1;
			// swap if necessary
			if (event.back_strip[small] > event.back_strip[big]) {
				small = 1;
				big = 0;
			}
			// small strip energy
			double se = event.back_energy[small];
			// big strip energy
			double be = event.back_energy[big];
			// small strip time
			double st = event.back_time[small] - ref_time;
			// big strip time
			double bt = event.back_time[big] - ref_time;
			// small strip
			unsigned short ss = event.back_strip[small] + dssd->Strip(0);
			// big strip
			unsigned short bs = event.back_strip[big] + dssd->Strip(0);
			// small strip module
			size_t smod = ss == dssd->Strip(0)+dssd->Strip(1)-1 ?
				ss / 16
				: (ss+1) / 16;
			// big strip module
			size_t bmod = bs == dssd->Strip(0) ?
				bs / 16
				: (bs-1) / 16;

			// fill single strip histograms
			time_energy_small[ss].Fill(se, st);
			time_small[ss].Fill(st);
			time_energy_big[bs].Fill(be, bt);
			time_big[bs].Fill(bt);
			// fill single module histograms
			time_energy_small_module[smod].Fill(se, st);
			time_small_module[smod].Fill(st);
			time_energy_big_module[bmod].Fill(be, bt);
			time_big_module[bmod].Fill(bt);
			// check cfd flag
			if ((event.cfd_flag & (1 << (small + 8))) == 0) {
				// fill single strip histograms
				time_energy_small_cfd[ss].Fill(se, st);
				time_small_cfd[ss].Fill(st);
				// fill single module histograms
				time_energy_small_module_cfd[smod].Fill(se, st);
				time_small_module_cfd[smod].Fill(st);
			} else {
				// fill single strip histograms
				time_energy_small_le[ss].Fill(se, st);
				time_small_le[ss].Fill(st);
				// fill single module histograms
				time_energy_small_module_le[smod].Fill(se, st);
				time_small_module_le[smod].Fill(st);
			}
			if ((event.cfd_flag & (1 << (big + 8))) == 0) {
				// fill single strip histograms
				time_energy_big_cfd[bs].Fill(be, bt);
				time_big_cfd[bs].Fill(bt);
				// fill single module histograms
				time_energy_big_module_cfd[bmod].Fill(be, bt);
				time_big_module_cfd[bmod].Fill(bt);
			} else {
				// fill single strip histograms
				time_energy_big_le[bs].Fill(be, bt);
				time_big_le[bs].Fill(bt);
				// fill single module histograms
				time_energy_big_module_le[bmod].Fill(be, bt);
				time_big_module_le[bmod].Fill(bt);
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save time-energy histograms of single strip and single hit
	for (TH2F &hist : time_energy_single) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_single_cfd) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_single_le) {
		hist.Write();
	}
	for (TH1F &hist : time_single) {
		hist.Write();
	}
	for (TH1F &hist : time_single_cfd) {
		hist.Write();
	}
	for (TH1F &hist : time_single_le) {
		hist.Write();
	}
	// save time-energy histograms of single strip and double hit small
	for (TH2F &hist : time_energy_small) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_small_cfd) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_small_le) {
		hist.Write();
	}
	for (TH1F &hist : time_small) {
		hist.Write();
	}
	for (TH1F &hist : time_small_cfd) {
		hist.Write();
	}
	for (TH1F &hist : time_small_le) {
		hist.Write();
	}
	// save time-energy histograms of single strip and double hit big
	for (TH2F &hist : time_energy_big) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_big_cfd) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_big_le) {
		hist.Write();
	}
	for (TH1F &hist : time_big) {
		hist.Write();
	}
	for (TH1F &hist : time_big_cfd) {
		hist.Write();
	}
	for (TH1F &hist : time_big_le) {
		hist.Write();
	}

	// save time-energy histograms of single module and single hit
	for (TH2F &hist : time_energy_single_module) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_single_module_cfd) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_single_module_le) {
		hist.Write();
	}
	for (TH1F &hist : time_single_module) {
		hist.Write();
	}
	for (TH1F &hist : time_single_module_cfd) {
		hist.Write();
	}
	for (TH1F &hist : time_single_module_le) {
		hist.Write();
	}
	// save time-energy histograms of single module and double hit small
	for (TH2F &hist : time_energy_small_module) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_small_module_cfd) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_small_module_le) {
		hist.Write();
	}
	for (TH1F &hist : time_small_module) {
		hist.Write();
	}
	for (TH1F &hist : time_small_module_cfd) {
		hist.Write();
	}
	for (TH1F &hist : time_small_module_le) {
		hist.Write();
	}
	// save time-energy histograms of single module and double hit big
	for (TH2F &hist : time_energy_big_module) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_big_module_cfd) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_big_module_le) {
		hist.Write();
	}
	for (TH1F &hist : time_big_module) {
		hist.Write();
	}
	for (TH1F &hist : time_big_module_cfd) {
		hist.Write();
	}
	for (TH1F &hist : time_big_module_le) {
		hist.Write();
	}
	// close files
	opf.Close();
	return 0;
}