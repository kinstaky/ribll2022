#include <iostream>

#include <TH1F.h>
#include <TRandom3.h>
#include <TString.h>
#include <TFile.h>
#include <TF1.h>

#include "include/defs.h"

using namespace ribll;

const int rate = 10'000;
const double mean_interval = 1e9 / rate;
const double dead_time = 8e3;

int main() {
	// output file name
	TString output_file_name = TString::Format(
		"%s%slive_time.root",
		kGenerateDataPath,
		kSimulateDir
	);
	// output file
	TFile output_file(output_file_name, "recreate");
	// histogram of interval distribution
	TH1F hist_interval("hi", "interval distribution", 1'000, 0, 100'000);

	TRandom3 generator(0);

	const int count = 30'000'000;
	const int count100 = count / 100 + 1;
	// show start
	printf("Simulating   0%%");
	fflush(stdout);
	double last_interval = 0.0;
	for (int i = 0; i < count; ++i) {
		// show process
		if (i % count100 == 0) {
			printf("\b\b\b\b%3d%%", i / count100);
			fflush(stdout);
		}
		double interval = mean_interval * -log(generator.Rndm());
		interval += last_interval;
		if (interval < dead_time) {
			last_interval = interval;
		} else {
			hist_interval.Fill(interval);
			last_interval = 0.0;
		}
	}
	printf("\b\b\b\b100%%\n");

	// fit
	TF1 fit_expo = TF1("f1", "expo", 10'000, 90'000);
	hist_interval.Fit(&fit_expo, "R+");
	std::cout << "Rate " << -1e9 * fit_expo.GetParameter(1) << "\n";

	// save histogram
	hist_interval.Write();
	// close file
	output_file.Close();
	return 0;
}