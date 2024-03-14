#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "include/event/dssd_event.h"
#include "include/detectors.h"

using namespace ribll;


void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] run d1_diff d2_diff\n"
		"  run               Set run number.\n"
		"  d1_diff           T0D1 energy difference.\n"
		"  d2_diff           T0D2 energy difference.\n"
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


bool SearchD2Front2Back1(
	const DssdFundamentalEvent &event,
	double diff,
	DssdMergeEvent &merge
) {
	// result
	bool found = false;
	// variables for convenient
	const int &fhit = event.front_hit;
	const int &bhit = event.back_hit;
	const unsigned short *fs = event.front_strip;
	const unsigned short *bs = event.back_strip;
	const double *fe = event.front_energy;
	const double *be = event.back_energy;
	const double *ft = event.front_time;

	// search for f2b1 event
	for (int i = 0; i < bhit; ++i) {
		for (int j = 0; j < fhit; ++j) {
			for (int k = j+1; k < fhit; ++k) {
				if (abs(fs[j]-fs[k]) == 1) continue;
				// energy difference
				double de = be[i] - fe[j] - fe[k];
				if (fabs(de) >= diff) continue;
				merge.hit = 2;
				merge.flag[0] = (0x100<<i) | (0x1<<j);
				merge.merge_tag[0] = 5;
				merge.energy[0] = fe[j];
				merge.time[0] = ft[j];
				merge.x[0] = fs[j];
				merge.y[0] = bs[i];
				merge.z[0] = 0.0;
				merge.flag[1] = (0x100<<i) | (0x1<<k);
				merge.merge_tag[1] = 5;
				merge.energy[1] = fe[k];
				merge.time[1] = ft[k];
				merge.x[1] = fs[k];
				merge.y[1] = bs[i];
				merge.z[1] = 0.0;
				found = true;
				break;
			}
			if (found) break;
		}
		if (found) break;
	}

	return found;
}


bool SearchD2Front1Back2(
	const DssdFundamentalEvent &event,
	double diff,
	DssdMergeEvent &merge
) {
	// result
	bool found = false;
	// variables for convenient
	const int &fhit = event.front_hit;
	const int &bhit = event.back_hit;
	const unsigned short *fs = event.front_strip;
	const unsigned short *bs = event.back_strip;
	const double *fe = event.front_energy;
	const double *be = event.back_energy;
	const double *ft = event.front_time;

	// search for f1b2 event
	for (int i = 0; i < fhit; ++i) {
		for (int j = 0; j < bhit; ++j) {
			for (int k = j+1; k < bhit; ++k) {
				if (abs(bs[j]-bs[k]) == 1) continue;
				// energy difference
				double de = fe[i] - be[j] - be[k];
				if (fabs(de) >= diff) continue;
				merge.hit = 2;
				merge.flag[0] = (0x1<<i) | (0x100<<j);
				merge.merge_tag[0] = 4;
				merge.energy[0] = fe[i] * be[j] / (be[j] + be[k]);
				merge.time[0] = ft[i];
				merge.x[0] = fs[i];
				merge.y[0] = bs[j];
				merge.z[0] = 0.0;
				merge.flag[1] = (0x1<<i) | (0x100<<k);
				merge.merge_tag[1] = 4;
				merge.energy[1] = fe[i] * be[k] / (be[j] + be[k]);
				merge.time[1] = ft[i];
				merge.x[1] = fs[i];
				merge.y[1] = bs[k];
				merge.z[1] = 0.0;
				found = true;
				break;
			}
			if (found) break;
		}
		if (found) break;
	}

	return found;
}


void CalculatePosition(
	Dssd *dssd,
	DssdMergeEvent &merge
) {
	for (int i = 0; i < 2; ++i) {
		auto position = dssd->CalculatePosition(merge.x[i], merge.y[i]);
		merge.x[i] = position.X();
		merge.y[i] = position.Y();
		merge.z[i] = position.Z();
	}
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

	if (pos_start+2 >= argc) {
		// positional arguments less than 1
		std::cerr << "Error: Miss run argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}

	int run = atoi(argv[pos_start]);
	double diff1 = atof(argv[pos_start+1]);
	double diff2 = atof(argv[pos_start+2]);

	// input T0D1 file name
	TString t0d1_file_name = TString::Format(
		"%s%st0d1-result-%s%04d.root",
		kGenerateDataPath,
		kNormalizeDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// input file
	TFile ipf(t0d1_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< t0d1_file_name << " failed.\n";
		return -1;
	}

	// input T0D2 file name
	TString t0d2_file_name = TString::Format(
		"%s%st0d2-result-%s%04d.root",
		kGenerateDataPath,
		kNormalizeDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// add friend
	ipt->AddFriend("d2=tree", t0d2_file_name);

	// T0D1 fundamental event
	DssdFundamentalEvent d1_event;
	// T0D2 fundamental event
	DssdFundamentalEvent d2_event;
	// setup input branches
	d1_event.SetupInput(ipt);
	d2_event.SetupInput(ipt, "d2.");

	// output T0D1 file name
	TString d1_merge_file_name = TString::Format(
		"%s%st0d1-merge-%ss1-%04d.root",
		kGenerateDataPath,
		kMergeDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// output T0D1 merge file
	TFile d1_merge_file(d1_merge_file_name, "recreate");
	// T0D1 merge tree
	TTree d1_tree("tree", "Merged binding events");
	// output T0D1 merge event
	DssdMergeEvent d1_merge;
	// setup output branches
	d1_merge.SetupOutput(&d1_tree);

	// output T0D2 file name
	TString d2_merge_file_name = TString::Format(
		"%s%st0d2-merge-%ss1-%04d.root",
		kGenerateDataPath,
		kMergeDir,
		tag.empty() ? "" : (tag+"-").c_str(),
		run
	);
	// output T0D2 merge file
	TFile d2_merge_file(d2_merge_file_name, "recreate");
	// T0D2 merge tree
	TTree d2_tree("tree", "Merged binding events");
	// output T0D2 merge event
	DssdMergeEvent d2_merge;
	// setup output branches
	d2_merge.SetupOutput(&d2_tree);

	T0d1 t0d1_detector(run, tag);
	T0d2 t0d2_detector(run, tag);

	// statistics
	long long d1_f1b2_count = 0;
	long long d1_f2ab2_count = 0;
	long long d1_f2b1_count = 0;
	long long d1_f2b2a_count = 0;
	long long d2_f1b2_count = 0;
	long long d2_f2b1_count = 0;

	// variable for convenient
	const int &fhit = d1_event.front_hit;
	const int &bhit = d1_event.back_hit;
	const unsigned short *fs = d1_event.front_strip;
	const unsigned short *bs = d1_event.back_strip;
	const double *fe = d1_event.front_energy;
	const double *be = d1_event.back_energy;
	const double *ft = d1_event.front_time;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Merging binding events   0%%");
	fflush(stdout);
	// loop
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get data
		ipt->GetEntry(entry);

		// initialize
		d1_merge.hit = 0;
		d2_merge.hit = 0;

		// search for D1 f1b2 event
		bool found_d1_f1b2 = false;
		for (int i = 0; i < fhit; ++i) {
			if (fe[i] < 15000.0) continue;
			// search for back twin events in energy range
			for (int j = 0; j < bhit; ++j) {
				if (be[j] < 10000.0) continue;
				for (int k = j+1; k < bhit; ++k) {
					if (be[k] > 10000.0) continue;
					// energy difference
					double de = fe[i] - be[j] - be[k];
					if (fabs(de) >= diff1) continue;
					// found
					d1_merge.hit = 2;
					d1_merge.flag[0] = (0x1<<i) | (0x100<<j);
					d1_merge.merge_tag[0] = 4;
					d1_merge.energy[0] = fe[i] * be[j] / (be[j]+be[k]);
					d1_merge.time[0] = ft[i];
					d1_merge.x[0] = fs[i];
					d1_merge.y[0] = bs[j];
					d1_merge.z[0] = 0.0;
					d1_merge.flag[1] = (0x1<<i) | (0x100<<k);
					d1_merge.merge_tag[1] = 4;
					d1_merge.energy[1] = fe[i] * be[k] / (be[j]+be[k]);
					d1_merge.time[1] = ft[i];
					d1_merge.x[1] = fs[i];
					d1_merge.y[1] = bs[k];
					d1_merge.z[1] = 0.0;
					found_d1_f1b2 = true;
					break;
				}
				if (found_d1_f1b2) break;
			}
			if (found_d1_f1b2) break;
		}
		if (found_d1_f1b2) {
			bool found_d2_f2b1 = SearchD2Front2Back1(d2_event, diff2, d2_merge);
			if (found_d2_f2b1) {
				CalculatePosition(&t0d1_detector, d1_merge);
				CalculatePosition(&t0d2_detector, d2_merge);
				d1_tree.Fill();
				d2_tree.Fill();
				++d1_f1b2_count;
				++d2_f2b1_count;
				continue;
			}
		}


		// search for D1 f2ab2 event
		bool found_d1_f2ab2 = false;
		for (int i = 0; i < fhit; ++i) {
			if (fe[i] < 10000.0) break;
			for (int j = i+1; j < fhit; ++j) {
				if (fe[j] > 10000.0) continue;
				if (abs(fs[i]-fs[j]) != 1) continue;
				// search for back events
				for (int k = 0; k < bhit; ++k) {
					double de1 = fe[i] - be[k];
					if (fabs(de1) >= diff1) continue;
					for (int l = k+1; l < bhit; ++l) {
						double de2 = fe[j] - be[l];
						if (fabs(de2) >= diff1) continue;
						// found
						d1_merge.hit = 2;
						d1_merge.flag[0] = (0x1<<i) | (0x100<<k);
						d1_merge.merge_tag[0] = 0;
						d1_merge.energy[0] = fe[i];
						d1_merge.time[0] = ft[i];
						d1_merge.x[0] = fs[i];
						d1_merge.y[0] = bs[k];
						d1_merge.z[0] = 0.0;
						d1_merge.flag[1] = (0x1<<j) | (0x100<<l);
						d1_merge.merge_tag[1] = 0;
						d1_merge.energy[1] = fe[j];
						d1_merge.time[1] = ft[j];
						d1_merge.x[1] = fs[j];
						d1_merge.y[1] = bs[l];
						d1_merge.z[1] = 0.0;
						found_d1_f2ab2 = true;
						break;
					}
					if (found_d1_f2ab2) break;
				}
				if (found_d1_f2ab2) break;
			}
			if (found_d1_f2ab2) break;
		}

		if (found_d1_f2ab2) {
			bool found_d2_f2b1 = SearchD2Front2Back1(d2_event, diff2, d2_merge);
			if (found_d2_f2b1) {
				CalculatePosition(&t0d1_detector, d1_merge);
				CalculatePosition(&t0d2_detector, d2_merge);
				d1_tree.Fill();
				d2_tree.Fill();
				++d1_f2ab2_count;
				++d2_f2b1_count;
				continue;
			}
		}

		// search for D1 f2b1 event
		bool found_d1_f2b1 = false;
		for (int i = 0; i < bhit; ++i) {
			if (be[i] < 15000.0) continue;
			// search for front twin events in energy range
			for (int j = 0; j < fhit; ++j) {
				if (fe[j] < 10000.0) continue;
				for (int k = j+1; k < fhit; ++k) {
					if (fe[k] > 10000.0) continue;
					// energy difference
					double de = be[i] - fe[j] - fe[k];
					if (fabs(de) >= diff1) continue;
					// found
					d1_merge.hit = 2;
					d1_merge.flag[0] = (0x100<<i) | (0x1<<j);
					d1_merge.merge_tag[0] = 5;
					d1_merge.energy[0] = fe[j];
					d1_merge.time[0] = ft[j];
					d1_merge.x[0] = fs[j];
					d1_merge.y[0] = bs[i];
					d1_merge.z[0] = 0.0;
					d1_merge.flag[1] = (0x100<<i) | (0x1<<k);
					d1_merge.merge_tag[1] = 5;
					d1_merge.energy[1] = fe[k];
					d1_merge.time[1] = ft[k];
					d1_merge.x[1] = fs[k];
					d1_merge.y[1] = bs[i];
					d1_merge.z[1] = 0.0;
					found_d1_f2b1 = true;
					break;
				}
				if (found_d1_f2b1) break;
			}
			if (found_d1_f2b1) break;
		}
		if (found_d1_f2b1) {
			bool found_d2_f1b2 = SearchD2Front1Back2(d2_event, diff2, d2_merge);
			if (found_d2_f1b2) {
				CalculatePosition(&t0d1_detector, d1_merge);
				CalculatePosition(&t0d2_detector, d2_merge);
				d1_tree.Fill();
				d2_tree.Fill();
				++d1_f2b1_count;
				++d2_f1b2_count;
				continue;
			}
		}


		// search for D1 f2b2a event
		bool found_d1_f2b2a = false;
		for (int i = 0; i < bhit; ++i) {
			if (be[i] < 10000.0) break;
			for (int j = i+1; j < bhit; ++j) {
				if (be[j] > 10000.0) continue;
				if (abs(bs[i]-bs[j]) != 1) continue;
				// search for front events
				for (int k = 0; k < fhit; ++k) {
					double de1 = be[i] - fe[k];
					if (fabs(de1) >= diff1) continue;
					for (int l = k+1; l < fhit; ++l) {
						double de2 = be[j] - fe[l];
						if (fabs(de2) >= diff1) continue;
						// found
						d1_merge.hit = 2;
						d1_merge.flag[0] = (0x100<<i) | (0x1<<k);
						d1_merge.merge_tag[0] = 0;
						d1_merge.energy[0] = fe[k];
						d1_merge.time[0] = ft[k];
						d1_merge.x[0] = fs[k];
						d1_merge.y[0] = bs[i];
						d1_merge.z[0] = 0.0;
						d1_merge.flag[1] = (0x100<<j) | (0x1<<l);
						d1_merge.merge_tag[1] = 0;
						d1_merge.energy[1] = fe[l];
						d1_merge.time[1] = ft[l];
						d1_merge.x[1] = fs[l];
						d1_merge.y[1] = bs[j];
						d1_merge.z[1] = 0.0;
						found_d1_f2b2a = true;
						break;
					}
					if (found_d1_f2b2a) break;
				}
				if (found_d1_f2b2a) break;
			}
			if (found_d1_f2b2a) break;
		}

		if (found_d1_f2b2a) {
			bool found_d2_f1b2 = SearchD2Front1Back2(d2_event, diff2, d2_merge);
			if (found_d2_f1b2) {
				CalculatePosition(&t0d1_detector, d1_merge);
				CalculatePosition(&t0d2_detector, d2_merge);
				d1_tree.Fill();
				d2_tree.Fill();
				++d1_f2b2a_count;
				++d2_f1b2_count;
				continue;
			}
		}

		d1_merge.hit = 0;
		d2_merge.hit = 0;
		d1_tree.Fill();
		d2_tree.Fill();
	}


	// show finish
	printf("\b\b\b\b100%%\n");

	// write tree and close files
	d1_merge_file.cd();
	d1_tree.Write();
	d1_merge_file.Close();
	d2_merge_file.cd();
	d2_tree.Write();
	d2_merge_file.Close();
	// close files
	ipf.Close();

	// print statistics
	std::cout << "D1 f1b2 " << d1_f1b2_count << "\n"
		<< "D1 f2ab2 " << d1_f2ab2_count << "\n"
		<< "D1 f2b1 " << d1_f2b1_count << "\n"
		<< "D1 f2b2a " << d1_f2b2a_count << "\n"
		<< "D2 f1b2 " << d2_f1b2_count << "\n"
		<< "D2 f2b1 " << d2_f2b1_count << "\n";
	return 0;
}