#include <iostream>

#include <TFile.h>
#include <TString.h>

#include "include/alignment.h"

using namespace ribll;

const int show_groups[] = {
	5, 10, 20, 50, 80, 100, 200, 500, 1000
};

int main(int argc, char **argv) {
	if (argc != 2) {
		std::cout << "Usage: " << argv[0] << " run" << std::endl;
		return -1;
	}

	// run number
	int run = atoi(argv[1]);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%salign-group-%04d.root",
		kGenerateDataPath, kShowDir, run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// match_rate-group graph
	TGraph rate_vs_group;
	// oversize-group graph
	TGraph oversize_vs_group;

	for (int group : show_groups) {
		Alignment align(
			run, group, 10'000'000,
			-10'000'000'000, 10'000'000'000
		);
		align.SetVerbose(false);
		if (align.Align()) {
			std::cerr << "Error: Align run " << run << " with groups "
				<< group << " failed.\n";
			opf.Close();
			return -1;
		}
		AlignStatistics statistics = align.GetStatistics();
		double match_rate =
			double(statistics.align_events) / double(statistics.VmeEvents());
		rate_vs_group.AddPoint(group, match_rate);
		oversize_vs_group.AddPoint(group, statistics.oversize_events);
	}
	// save graphs
	opf.cd();
	rate_vs_group.Write("grate");
	oversize_vs_group.Write("gover");
	// close file
	opf.Close();
	return 0;
}