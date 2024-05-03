#include <iostream>

#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TString.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/channel_event.h"

using namespace ribll;

class Spectrum {
public:

	Spectrum(int n): n_(n) {}

	double operator()(double *x, double *par) {
		double result = 0.0;
		result += par[0] + par[1] * x[0];
		for (int i = 0; i < n_; ++i) {
			result += par[3*i+2] * exp(-0.5 * pow((x[0] - par[3*i+3]) / par[3*i+4], 2.0));
		}
		return result;
	}

private:
	int n_;
};


int main() {
	TChain chain("tree", "chain");
	for (int run = 618; run <= 716; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		chain.AddFile(TString::Format(
			"%s%sBe8-4He-4He-%04d.root/tree",
			kGenerateDataPath,
			kChannelDir,
			run
		));
    }
	// input data
	ChannelEvent event;
	// setup input branches
	event.SetupInput(&chain);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s2a.root",
		kGenerateDataPath,
		kSpectrumDir
	);
	// output file
	TFile output_file(output_file_name, "recreate");

	// histogram of Q value
	TH1F hist("h", "2 alpha", 1000, -1, 4);
	hist.SetLineColor(kBlack);

	for (long long entry = 0; entry < chain.GetEntriesFast(); ++entry) {
		chain.GetEntry(entry);

		double ex = event.daughter_energy[0] + event.daughter_energy[1]
			- event.parent_energy - 0.092;
		hist.Fill(ex);
	}

	// fit total Q spectrum
	Spectrum spectrum(1);
	// fit histograms
	TF1 fit_spectrum("fs", spectrum, -0.1, 0.08, 5);
	double fit_initial_parameters[5] = {
		0.0, 0.5,
		1000.0, 0.0, 0.05
	};
	fit_spectrum.SetParameters(fit_initial_parameters);
	// fit_spectrum.SetParLimits(0, 0.1, 100.0);
	// fit_spectrum.SetParLimits(1, -13.0, -11.0);
	fit_spectrum.SetParLimits(2, 0.0, 10000.0);
	fit_spectrum.SetParLimits(3, -0.5, 0.5);
	// fit_spectrum.SetParLimits(4, 0.0, 0.5);
	fit_spectrum.SetNpx(1000);
	hist.Fit(&fit_spectrum, "R+");

	double parameters[5];
	fit_spectrum.GetParameters(parameters);
	TF1 *fit_specs[1];
	for (int i = 0; i < 1; ++i) {
		fit_specs[i] = new TF1(
			TString::Format("f%d", i), "gaus", -0.1, 0.1
		);
		fit_specs[i]->SetParameters(parameters+3*i+2);
		fit_specs[i]->SetLineColor(kBlue);
		fit_specs[i]->SetNpx(200);
		hist.GetListOfFunctions()->Add(fit_specs[i]);
	}
	// TF1 *fit_bg = new TF1("fbg", "pol1", 7, 11);
	// fit_bg->SetParameters(parameters);
	// fit_bg->SetLineColor(kBlue);
	// hist.GetListOfFunctions()->Add(fit_bg);

	// save
	hist.Write();
	// close files
	output_file.Close();
	return 0;
}