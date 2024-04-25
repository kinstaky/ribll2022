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
		if (x[0] > 8.0) result += par[0] + par[1] * x[0];
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
			"%s%sC12-4He-4He-4He-%04d.root/tree",
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
		"%s%s3a.root",
		kGenerateDataPath,
		kSpectrumDir
	);
	// output file
	TFile output_file(output_file_name, "recreate");

	// histogram of Q value
	TH1F hist("h", "3 alpha", 200, 7, 37);
	hist.SetLineColor(kBlack);

	for (long long entry = 0; entry < chain.GetEntriesFast(); ++entry) {
		chain.GetEntry(entry);

		double ex = event.daughter_energy[0] + event.daughter_energy[1]
			+ event.daughter_energy[2] - event.parent_energy + 7.275;
		hist.Fill(ex);
	}

	// fit total Q spectrum
	Spectrum spectrum(2);
	// fit histograms
	TF1 fit_spectrum("fs", spectrum, 7, 11, 8);
	double fit_initial_parameters[8] = {
		0.0, 0.5,
		100.0, 7.6, 0.05,
		100.0, 9.6, 0.05
	};
	fit_spectrum.SetParameters(fit_initial_parameters);
	// fit_spectrum.SetParLimits(0, 0.1, 100.0);
	// fit_spectrum.SetParLimits(1, -13.0, -11.0);
	fit_spectrum.SetParLimits(2, 0.0, 100.0);
	fit_spectrum.SetParLimits(3, 7.0, 8.0);
	fit_spectrum.SetParLimits(4, 0.0, 0.5);
	fit_spectrum.SetParLimits(5, 0.1, 200.0);
	fit_spectrum.SetParLimits(6, 9.0, 10.0);
	fit_spectrum.SetParLimits(7, 0.0, 0.5);
	fit_spectrum.SetNpx(1000);
	hist.Fit(&fit_spectrum, "R+");

	double parameters[8];
	fit_spectrum.GetParameters(parameters);
	TF1 *fit_specs[2];
	for (int i = 0; i < 2; ++i) {
		fit_specs[i] = new TF1(
			TString::Format("f%d", i), "gaus", 7, 17
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