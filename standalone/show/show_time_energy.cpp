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
		"%s%s%s-time-energy-%s%s%04u-%04u.root",
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
	// time-energy 2D histogram of single strip
	std::vector<TH2F> time_energy;
	for (size_t side = 0;  side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_energy.emplace_back(
				TString::Format("%cte%ld", "fb"[side], i),
				"time energy",
				1000, 0, 60000,
				2000, normalize ? -500 : -200,  normalize ? 500 : 800
			);
		}
	}
	// time-energy 2D histogram of single strip with valid CFD
	std::vector<TH2F> time_energy_cfd;
	for (size_t side = 0;  side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_energy_cfd.emplace_back(
				TString::Format("%ctec%ld", "fb"[side], i),
				"time energy with CFD",
				1000, 0, 60000,
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time-energy 2D histogram of single strip with leading edge detected
	std::vector<TH2F> time_energy_le;
	for (size_t side = 0;  side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			time_energy_le.emplace_back(
				TString::Format("%ctel%ld", "fb"[side], i),
				"time energy with LE",
				1000, 0, 60000,
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time 1D histogram of single strip
	std::vector<TH1F> hist_time;
	for (size_t side = 0; side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			hist_time.emplace_back(
				TString::Format("%ct%ld", "fb"[side], i),
				"time",
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time 1D histogram of single strip with CFD
	std::vector<TH1F> hist_time_cfd;
	for (size_t side = 0; side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			hist_time_cfd.emplace_back(
				TString::Format("%ctc%ld", "fb"[side], i),
				"time with CFD",
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time 1D histogram of single strip with LE
	std::vector<TH1F> hist_time_le;
	for (size_t side = 0; side < 2; ++side) {
		for (size_t i = 0; i < dssd->Strip(side); ++i) {
			hist_time_le.emplace_back(
				TString::Format("%ctl%ld", "fb"[side], i),
				"time with LE",
				2000, normalize ? -500 : -200, normalize ? 500 : 800
			);
		}
	}
	// time-energy 2D histogram of single module
	std::vector<TH2F> time_energy_module;
	for (size_t i = 0; i < 8; ++i) {
		time_energy_module.emplace_back(
			TString::Format("mte%ld", i),
			"time-energy of single module",
			1000, 0, 60000,
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time-energy 2D histogram of single module with CFD
	std::vector<TH2F> time_energy_module_cfd;
	for (size_t i = 0; i < 8; ++i) {
		time_energy_module_cfd.emplace_back(
			TString::Format("mtec%ld", i),
			"time-energy of single module with CFD",
			1000, 0, 60000,
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time-energy 2D histogram of single module with LE
	std::vector<TH2F> time_energy_module_le;
	for (size_t i = 0; i < 8; ++i) {
		time_energy_module_le.emplace_back(
			TString::Format("mtel%ld", i),
			"time-energy of single module with LE",
			1000, 0, 60000,
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time 1D histogram of single module
	std::vector<TH1F> time_module;
	for (size_t i = 0; i < 8; ++i) {
		time_module.emplace_back(
			TString::Format("mt%ld", i),
			"time of single module",
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time 1D histogram of single module with CFD
	std::vector<TH1F> time_module_cfd;
	for (size_t i = 0; i < 8; ++i) {
		time_module_cfd.emplace_back(
			TString::Format("mtc%ld", i),
			"time of single module with CFD",
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time 1D histogram of single module with LE
	std::vector<TH1F> time_module_le;
	for (size_t i = 0; i < 8; ++i) {
		time_module_le.emplace_back(
			TString::Format("mtl%ld", i),
			"time of single module with LE",
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time-energy 2D histogram of single side
	std::vector<TH2F> time_energy_side;
	for (size_t i = 0; i < 2; ++i) {
		time_energy_side.emplace_back(
			TString::Format("ste%ld", i),
			"time-energy of single side",
			1000, 0, 60000,
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time-eneryg 2D histogram of single side with CFD
	std::vector<TH2F> time_energy_side_cfd;
	for (size_t i = 0; i < 2; ++i) {
		time_energy_side_cfd.emplace_back(
			TString::Format("stec%ld", i),
			"time-energy of single side with CFD",
			1000, 0, 60000,
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time-energy 2D histogram of single side with LE
	std::vector<TH2F> time_energy_side_le;
	for (size_t i = 0; i < 2; ++i) {
		time_energy_side_le.emplace_back(
			TString::Format("stel%ld", i),
			"time-energy of single side with LE",
			1000, 0, 60000,
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time 1D histogram of single side
	std::vector<TH1F> time_side;
	for (size_t i = 0; i < 2; ++i) {
		time_side.emplace_back(
			TString::Format("st%ld", i),
			"time of single side",
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time 1D histogram of single side with CFD
	std::vector<TH1F> time_side_cfd;
	for (size_t i = 0; i < 2; ++i) {
		time_side_cfd.emplace_back(
			TString::Format("stc%ld", i),
			"time of single side with CFD",
			2000, normalize ? -500 : -200, normalize ? 500 : 800
		);
	}
	// time 1D histgoram of single side with LE
	std::vector<TH1F> time_side_le;
	for (size_t i = 0; i < 2; ++i) {
		time_side_le.emplace_back(
			TString::Format("stl%ld", i),
			"time of single side with LE",
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
		// fill front events
		for (unsigned short i = 0; i < event.front_hit; ++i) {
			// fill single strip histograms
			time_energy[event.front_strip[i]].Fill(
				event.front_energy[i],
				event.front_time[i] - ref_time
			);
			hist_time[event.front_strip[i]].Fill(
				event.front_time[i] - ref_time
			);
			// fill single module histograms
			time_energy_module[event.front_strip[i]/16].Fill(
				event.front_energy[i],
				event.front_time[i] - ref_time
			);
			time_module[event.front_strip[i]/16].Fill(
				event.front_time[i] - ref_time
			);
			// fill single side histograms
			time_energy_side[0].Fill(
				event.front_energy[i],
				event.front_time[i] - ref_time
			);
			time_side[0].Fill(
				event.front_time[i] - ref_time
			);
			if ((event.cfd_flag & (1 << i)) == 0) {
				// fill single strip histograms
				time_energy_cfd[event.front_strip[i]].Fill(
					event.front_energy[i],
					event.front_time[i] - ref_time
				);
				hist_time_cfd[event.front_strip[i]].Fill(
					event.front_time[i] - ref_time
				);
				// fill single module histograms
				time_energy_module_cfd[event.front_strip[i]/16].Fill(
					event.front_energy[i],
					event.front_time[i] - ref_time
				);
				time_module_cfd[event.front_strip[i]/16].Fill(
					event.front_time[i] - ref_time
				);
				// fill single side histograms
				time_energy_side_cfd[0].Fill(
					event.front_energy[i],
					event.front_time[i] - ref_time
				);
				time_side_cfd[0].Fill(
					event.front_time[i] - ref_time
				);
			} else {
				// fill single strip histograms
				time_energy_le[event.front_strip[i]].Fill(
					event.front_energy[i],
					event.front_time[i] - ref_time
				);
				hist_time_le[event.front_strip[i]].Fill(
					event.front_time[i] - ref_time
				);
				// fill single module histograms
				time_energy_module_le[event.front_strip[i]/16].Fill(
					event.front_energy[i],
					event.front_time[i] - ref_time
				);
				time_module_le[event.front_strip[i]/16].Fill(
					event.front_time[i] - ref_time
				);
				// fill single side histograms
				time_energy_side_le[0].Fill(
					event.front_energy[i],
					event.front_time[i] - ref_time
				);
				time_side_le[0].Fill(
					event.front_time[i] - ref_time
				);
			}
		}
		// fill back events
		for (unsigned short i = 0; i < event.back_hit; ++i) {
			size_t strip = dssd->Strip(0) + event.back_strip[i];
			// fill single strip histograms
			time_energy[strip].Fill(
				event.back_energy[i],
				event.back_time[i] - ref_time
			);
			hist_time[strip].Fill(
				event.back_time[i] - ref_time
			);
			// fill single module histograms
			time_energy_module[strip/16].Fill(
				event.back_energy[i],
				event.back_time[i] - ref_time
			);
			time_module[strip/16].Fill(
				event.back_time[i] - ref_time
			);
			// fill single side histograms
			time_energy_side[1].Fill(
				event.back_energy[i],
				event.back_time[i] - ref_time
			);
			time_side[1].Fill(
				event.back_time[i] - ref_time
			);
			if (((event.cfd_flag & (1 << (i + 8)))) == 0) {
				// fill single strip histograms
				time_energy_cfd[strip].Fill(
					event.back_energy[i],
					event.back_time[i] - ref_time
				);
				hist_time_cfd[strip].Fill(
					event.back_time[i] - ref_time
				);
				// fill single module histograms
				time_energy_module_cfd[strip/16].Fill(
					event.back_energy[i],
					event.back_time[i] - ref_time
				);
				time_module_cfd[strip/16].Fill(
					event.back_time[i] - ref_time
				);
				// fill single side histograms
				time_energy_side_cfd[1].Fill(
					event.back_energy[i],
					event.back_time[i] - ref_time
				);
				time_side_cfd[1].Fill(
					event.back_time[i] - ref_time
				);
			} else {
				// fill single strip histograms
				time_energy_le[strip].Fill(
					event.back_energy[i],
					event.back_time[i] - ref_time
				);
				hist_time_le[strip].Fill(
					event.back_time[i] - ref_time
				);
				// fill single module histograms
				time_energy_module_le[strip/16].Fill(
					event.back_energy[i],
					event.back_time[i] - ref_time
				);
				time_module_le[strip/16].Fill(
					event.back_time[i] - ref_time
				);
				// fill single side histograms
				time_energy_side_le[1].Fill(
					event.back_energy[i],
					event.back_time[i] - ref_time
				);
				time_side_le[1].Fill(
					event.back_time[i] - ref_time
				);
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save time-energy histograms of single strip
	for (TH2F &hist : time_energy) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_cfd) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_le) {
		hist.Write();
	}
	for (TH1F &hist : hist_time) {
		hist.Write();
	}
	for (TH1F &hist : hist_time_cfd) {
		hist.Write();
	}
	for (TH1F &hist : hist_time_le) {
		hist.Write();
	}
	// save time-energy histograms of single module
	for (TH2F &hist : time_energy_module) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_module_cfd) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_module_le) {
		hist.Write();
	}
	for (TH1F &hist : time_module) {
		hist.Write();
	}
	for (TH1F &hist : time_module_cfd) {
		hist.Write();
	}
	for (TH1F &hist : time_module_le) {
		hist.Write();
	}
	// save time-energy histograms of single side
	for (TH2F &hist : time_energy_side) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_side_cfd) {
		hist.Write();
	}
	for (TH2F &hist : time_energy_side_le) {
		hist.Write();
	}
	for (TH1F &hist : time_side) {
		hist.Write();
	}
	for (TH1F &hist : time_side_cfd) {
		hist.Write();
	}
	for (TH1F &hist : time_side_le) {
		hist.Write();
	}
	// close files
	opf.Close();
	return 0;
}