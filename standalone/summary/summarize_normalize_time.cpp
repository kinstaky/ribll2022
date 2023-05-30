#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TString.h>

#include "include/detectors.h"

using namespace ribll;

int main(int argc, char **argv) {
	if (argc < 4) {
		std::cout << "Usage: " << argv[0] << " run end_run detector\n"
			"  run               start run, inclusive.\n"
			"  end_run           end run, inclusive.\n"
			"  detector          detector name.\n";
		return -1;
	}
	// start run
	unsigned int run  = atoi(argv[1]);
	// end run, inclusive
	unsigned int end_run = atoi(argv[2]);
	// detector name
	std::string detector_name(argv[3]);
	// pointer to DSSD
	std::shared_ptr<Dssd> dssd = CreateDssd(detector_name, run, "");

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-norm-time-%04u-%04u.root",
		kGenerateDataPath,
		kSummaryDir,
		detector_name.c_str(),
		run,
		end_run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// parameters of front side
	std::vector<TGraph> front_params;
	for (size_t i = 0; i < dssd->Strip(0); ++i) {
		front_params.emplace_back();
	}
	// parameters of back side
	std::vector<TGraph> back_params;
	for (size_t i = 0; i < dssd->Strip(1); ++i) {
		back_params.emplace_back();
	}

	for (unsigned int i = run; i <= end_run; ++i) {
		if (i == 628) continue;
		std::ifstream fin(TString::Format(
			"%s%s%s-norm-param-%04u.txt",
			kGenerateDataPath,
			kTimeDir,
			detector_name.c_str(),
			i
		).Data());
		if (!fin.good()) {
			std::cout << "Warning: Read parameters of run "
				<< i << " failed.\n";
			continue;
		}
		for (size_t j = 0; j < dssd->Strip(0); ++j) {
			double param;
			fin >> param;
			front_params[j].AddPoint(i, param);
		}
		for (size_t j = 0; j < dssd->Strip(1); ++j) {
			double param;
			fin >> param;
			back_params[j].AddPoint(i, param);
		}
		// close file
		fin.close();
	}
	// save graphs
	for (size_t i = 0; i < dssd->Strip(0); ++i) {
		front_params[i].Write(TString::Format("f%ld", i));
	}
	for (size_t i = 0; i < dssd->Strip(1); ++i) {
		back_params[i].Write(TString::Format("b%ld", i));
	}
	// close files
	opf.Close();
	return 0;
}