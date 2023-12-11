#include "include/mapper/trace_mapper.h"

#include <iostream>
#include <vector>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include "include/defs.h"
#include "include/event/dssd_event.h"

namespace ribll {


TraceMapper::TraceMapper(unsigned int run, unsigned int crate)
: run_(run)
, crate_(crate) {

}

int TraceMapper::Map() {
	if (crate_ == 0) {
		// map tafd trace
		for (int i = 2; i < 4; ++i) {
			// list of decode entries
			std::vector<long long> input_entries;
			// input entry file name
			TString entry_file_name;
			entry_file_name.Form(
				"%s%stafd%d-merge-ta-%04u.root",
				kGenerateDataPath,
				kMergeDir,
				i,
				run_
			);
			// input entry file
			TFile entry_file(entry_file_name, "read");
			// input entry tree
			TTree *entry_tree = (TTree*)entry_file.Get("tree");
			if (!entry_tree) {
				std::cerr << "Error: Get tree from "
					<< entry_file_name << " failed.\n";
				return -1;
			}
			// input merge event
			AdssdMergeEvent event;
			// setup input entry tree branches
			event.SetupInput(entry_tree);
			// total number of entries
			long long entries = entry_tree->GetEntries();
			// read entries
			for (long long entry = 0; entry < entries; ++entry) {
				entry_tree->GetEntry(entry);
				if (event.hit == 1) {
					input_entries.push_back(event.decode_entry[0]);
				} else {
					input_entries.push_back(-1);
				}
			}
			// close file
			entry_file.Close();

			// input decode file name
			TString decode_file_name;
			decode_file_name.Form(
				"%s%s_R%04u.root",
				kCrate0Path, kCrate0FileName, run_
			);
			// input decode file
			TFile ipf(decode_file_name, "read");
			// input decode tree
			TTree *ipt = (TTree*)ipf.Get("tree");
			// input trace points
			unsigned short ltra;
			// input trace data
			unsigned short trace[16384];
			// setup input branches
			ipt->SetBranchAddress("ltra", &ltra);
			ipt->SetBranchAddress("data", trace);

			// output trace file name
			TString trace_file_name;
			trace_file_name.Form(
				"%s%stafd%d-trace-ta-%04u.root",
				kGenerateDataPath,
				kTraceDir,
				i,
				run_
			);
			// output trace file
			TFile opf(trace_file_name, "recreate");
			// output trace tree
			TTree opt("tree", "trace");
			// output data
			int points;
			// setup output branches
			opt.Branch("point", &points, "p/I");
			opt.Branch("trace", trace, "t[p]/s");

			// 1/100 of entries, for showing process
			long long entry100 = entries / 100;
			// show start
			printf("Mapping trace of tafd%d   0%%", i);
			fflush(stdout);
			for (long long entry = 0; entry < entries; ++entry) {
				// show process
				if (entry % entry100 == 0) {
					printf("\b\b\b\b%3lld%%", entry / entry100);
					fflush(stdout);
				}
				if (input_entries[entry] >= 0) {
					ipt->GetEntry(input_entries[entry]);
					points = ltra;
				} else {
					points = 0;
				}
				opt.Fill();
			}
			// show finish
			printf("\b\b\b\b100%%\n");
			// save tree
			opt.Write();
			// close files
			opf.Close();
			ipf.Close();
		}

		// // map csi trace
		// // list of decode entries
		// for (int i = 0; i < 6; ++i) {
		// 	std::vector<long long> input_entries;
		// 	// input entry file name
		// 	TString entry_file_name;
		// 	entry_file_name.Form(
		// 		"%s%staf%d-telescope-ta-%04u.root",
		// 		kGenerateDataPath,
		// 		kTelescopeDir,
		// 		i,
		// 		run_
		// 	);
		// 	// input entry file
		// 	TFile entry_file(entry_file_name, "read");
		// 	// input entry tree
		// 	TTree *entry_tree = (TTree*)entry_file.Get("tree");
		// 	if (!entry_tree) {
		// 		std::cerr << "Error: Get tree from "
		// 			<< entry_file_name << " failed.\n";
		// 		return -1;
		// 	}
		// 	// input event number
		// 	unsigned short num;
		// 	// input csi decode entry
		// 	long long csi_decode_entry[4];
		// 	// setup input entry tree branches
		// 	entry_tree->SetBranchAddress("num", &num);
		// 	entry_tree->SetBranchAddress("csi_decode_entry", csi_decode_entry);
		// 	// total number of entries
		// 	long long entries = entry_tree->GetEntries();
		// 	// read entries
		// 	for (long long entry = 0; entry < entries; ++entry) {
		// 		entry_tree->GetEntry(entry);
		// 		if (num == 1) {
		// 			input_entries.push_back(csi_decode_entry[0]);
		// 		} else {
		// 			input_entries.push_back(-1);
		// 		}
		// 	}
		// 	// close file
		// 	entry_file.Close();

		// 	// input decode file name
		// 	TString decode_file_name;
		// 	decode_file_name.Form(
		// 		"%s%s_R%04u.root",
		// 		kCrate0Path, kCrate0FileName, run_
		// 	);
		// 	// input decode file
		// 	TFile ipf(decode_file_name, "read");
		// 	// input decode tree
		// 	TTree *ipt = (TTree*)ipf.Get("tree");
		// 	// input trace points
		// 	unsigned short ltra;
		// 	// input trace data
		// 	unsigned short trace[16384];
		// 	// setup input branches
		// 	ipt->SetBranchAddress("ltra", &ltra);
		// 	ipt->SetBranchAddress("data", trace);

		// 	// output trace file name
		// 	TString trace_file_name;
		// 	trace_file_name.Form(
		// 		"%s%stafcsi%d%d-trace-ta-%04u.root",
		// 		kGenerateDataPath,
		// 		kTraceDir,
		// 		i*2,
		// 		i*2+1,
		// 		run_
		// 	);
		// 	// output trace file
		// 	TFile opf(trace_file_name, "recreate");
		// 	// output trace tree
		// 	TTree opt("tree", "trace");
		// 	// setup output branches
		// 	opt.Branch("point", &ltra, "p/s");
		// 	opt.Branch("trace", trace, "t[p]/s");

		// 	// 1/100 of entries, for showing process
		// 	long long entry100 = entries / 100;
		// 	// show start
		// 	printf("Mapping trace of tafcsi%d%d   0%%", i*2, i*2+1);
		// 	fflush(stdout);
		// 	for (long long entry = 0; entry < entries; ++entry) {
		// 		// show process
		// 		if (entry % entry100 == 0) {
		// 			printf("\b\b\b\b%3lld%%", entry / entry100);
		// 			fflush(stdout);
		// 		}
		// 		if (input_entries[entry] >= 0) {
		// 			ipt->GetEntry(input_entries[entry]);
		// 		} else {
		// 			ltra = 0;
		// 		}
		// 		opt.Fill();
		// 	}
		// 	// show finish
		// 	printf("\b\b\b\b100%%\n");
		// 	// save tree
		// 	opt.Write();
		// 	// close files
		// 	opf.Close();
		// 	ipf.Close();
		// }
	} else {
		std::cerr << "Error: TraceMapper for crate " << crate_
			<< " is not implemented yet.\n";
		return -1;
	}
	return 0;
}

}	// namespace ribll