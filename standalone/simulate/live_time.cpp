#include <iostream>
#include <cmath>

#include <TH1F.h>
#include <TRandom3.h>
#include <TString.h>
#include <TFile.h>
#include <TF1.h>

#include "include/defs.h"

using namespace ribll;

const int rate = 10'000;
const double mean_interval = 1e9 / rate;
const double dead_time = 4e3;

int Factorial(int n) {
	if (n == 0) return 1;
	else if (n == 1) return 1;
	else return n * Factorial(n-1);
}

double FitExtended(double *x, double *par) {
	double result = 0.0;
	double t = x[0];
	double a = par[0];
	double r = par[1];
	int i = 1;
	while (t > i * dead_time) {
		double tmp = t - i * dead_time;
		tmp = pow(r, i) / Factorial(i-1) * pow(tmp, i-1) * exp(-i*r*dead_time);
		tmp = i % 2 == 0 ? -tmp : tmp;
		result += tmp;
		++i;
	}
	return a * result;
}

int main() {
	// output file name
	TString output_file_name = TString::Format(
		"%s%slive_time.root",
		kGenerateDataPath,
		kSimulateDir
	);
	// output file
	TFile output_file(output_file_name, "recreate");
	// histogram of non-extended model interval distribution
	TH1F hist_interval_non(
		"hin", "interval distribution of non-extended model",
		1'000, 0, 100'000
	);
	// histogram of extended interval distribution
	TH1F hist_interval_extended(
		"hie", "interval distribution of extended model",
		1'000, 0, 100'000
	);

	TRandom3 generator(0);

	const int count = 30'000'000;
	const int count100 = count / 100 + 1;
	// show start
	printf("Simulating   0%%");
	fflush(stdout);
	// last interval of non-extended model
	double last_interval_non = 0.0;
	// last interval of extended model
	double last_interval_extended = 0.0;
	for (int i = 0; i < count; ++i) {
		// show process
		if (i % count100 == 0) {
			printf("\b\b\b\b%3d%%", i / count100);
			fflush(stdout);
		}
		double interval = mean_interval * -log(generator.Rndm());
		double interval_non = interval + last_interval_non;
		if (interval_non < dead_time) {
			last_interval_non = interval_non;
		} else {
			hist_interval_non.Fill(interval_non);
			last_interval_non = 0.0;
		}
		double interval_extended = interval + last_interval_extended;
		if (interval < dead_time) {
			last_interval_extended = interval_extended;
		} else {
			hist_interval_extended.Fill(interval_extended);
			last_interval_extended = 0.0;
		}
	}
	printf("\b\b\b\b100%%\n");

	// fit non-extended model
	TF1 fit_non = TF1("fn", "expo", 10'000, 90'000);
	hist_interval_non.Fit(&fit_non, "R+");
	std::cout << "Non-extended rate "
		<< -1e9 * fit_non.GetParameter(1) << "\n";

	// fit extended model
	TF1 fit_extended = TF1("fe", FitExtended, 4'000, 100'000, 2);
	fit_extended.SetNpx(10000);
	fit_extended.SetParameter(0, 1e11);
	fit_extended.SetParameter(1, 1e-5);
	hist_interval_extended.Fit(&fit_extended, "R+");
	std::cout << "Extended rate "
		<< 1e9 * fit_extended.GetParameter(1) << "\n";

	// save histogram
	hist_interval_non.Write();
	hist_interval_extended.Write();
	// close file
	output_file.Close();
	return 0;
}