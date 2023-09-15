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

constexpr double be10_mass = 10.0113403769;
constexpr double he4_mass = 4.0015060943;
constexpr double h2_mass = 2.0135531980;
constexpr double c14_mass = 13.9999505089;

constexpr double u = 931.494;

double QSpectrum(double *x, double *par) {
	double result = 0.0;
	for (int i = 0; i < 4; ++i) {
		result += par[3*i] * exp(-0.5*pow((x[0]-par[3*i+1])/par[3*i+2], 2.0));
	}
	return result;
}

int main() {
	TChain input_chain("tree", "input channel chain");
	for (unsigned int run = 618; run <= 716; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		input_chain.AddFile(TString::Format(
			"%s%sC14-10Be-4He-2H-%04u.root/tree",
			kGenerateDataPath,
			kChannelDir,
			run
		));
	}
	ChannelEvent channel;
	// setup input tree
	channel.SetupInput(&input_chain);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%sthreebody.root",
		kGenerateDataPath,
		kSpectrumDir
	);
	// output file
	TFile output_file(output_file_name, "recreate");

	// histogram of Q value
	TH1F hist_q("hq", "Q value", 60, -23, -8);
	hist_q.SetLineColor(kBlack);
	// histogram of 14C decay to all state of 10Be
	TH1F c_spec_all("hsca", "spectrum of 14C", 100, 10, 40);
	// spectrum of 14C decay to 10Be ground state
	TH1F c_spec_0("hsc0", "spectrum of 14C to 10Be ground state", 50, 10, 40);
	// spectrum of 14C decay to 10Be 3.5 MeV state
	TH1F c_spec_1("hsc1", "spectrum of 14C to 10Be 3.5MeV state", 60, 10, 40);
	// spectrum of 14C decay to 10Be 6 MeV state
	TH1F c_spec_2("hsc2", "spectrum of 14C to 10Be 6MeV state", 60, 10, 40);
	// spectrum of 14C decay to 10Be 7.5 MeV state
	TH1F c_spec_3("hsc3", "spectrum of 14C to 10Be 7.5MeV state", 60, 10, 40);

	for (long long entry = 0; entry < input_chain.GetEntriesFast(); ++entry) {
		input_chain.GetEntry(entry);
		double q = channel.daughter_energy[0] + channel.daughter_energy[1]
			+ channel.recoil_energy - channel.parent_energy;
		hist_q.Fill(q);
		// 10Be excited energy
		double be_excited = -q - 11.3;
		// double be_excited = 0.0;
		// if (q > -13.5) be_excited = 0.0;
		// else if (q > -16.5) be_excited = 3.5;
		// else if (q > -19) be_excited = 6.0;
		// else be_excited = 7.5;
		// calculate 14C excited spectrum
		ROOT::Math::XYZVector be_momentum(
			channel.daughter_px[0],
			channel.daughter_py[0],
			channel.daughter_pz[0]
		);
		double be_energy = sqrt(
			pow(be_momentum.R(), 2.0)
			+ pow(be10_mass*u+be_excited, 2.0)
		);
		ROOT::Math::XYZVector he_momentum(
			channel.daughter_px[1],
			channel.daughter_py[1],
			channel.daughter_pz[1]
		);
		double he_energy = sqrt(
			pow(he_momentum.R(), 2.0)
			+ pow(he4_mass*u, 2.0)
		);
		ROOT::Math::XYZVector c_momentum = be_momentum + he_momentum;
		double c_energy = be_energy + he_energy;
		double c_excited_mass = sqrt(
			pow(c_energy, 2.0)
			- pow(c_momentum.R(), 2.0)
		);
		double c_excited = c_excited_mass - c14_mass*u;

		if (channel.parent_energy < 378 || channel.parent_energy > 390) continue;
		c_spec_all.Fill(c_excited);
		if (q < -10.5 && q > -13) c_spec_0.Fill(c_excited);
		if (q < -14 && q > -15.5) c_spec_1.Fill(c_excited);
		if (q < -16.3 && q > -18.5) c_spec_2.Fill(c_excited);
		if (q < -20.0 && q > -21.0) c_spec_3.Fill(c_excited);
	}

	// fit histograms
	TF1 fit_q_spectrum("fq", QSpectrum, -23, -8, 12);
	double fit_q_initial_parameters[12] = {
		20.0, -11.5, 1.0,
		60.0, -15.0, 1.0,
		50.0, -17.5, 1.0,
		10.0, -19.0, 0.5
	};
	fit_q_spectrum.SetParameters(fit_q_initial_parameters);
	// set limits to A
	fit_q_spectrum.SetParLimits(0, 0.1, 100.0);
	fit_q_spectrum.SetParLimits(3, 0.1, 100.0);
	fit_q_spectrum.SetParLimits(6, 0.1, 100.0);
	fit_q_spectrum.SetParLimits(9, 1.0, 100.0);
	// set limits to mean
	fit_q_spectrum.SetParLimits(1, -13.0, -11.0);
	fit_q_spectrum.SetParLimits(4, -16.0, -13.5);
	fit_q_spectrum.SetParLimits(7, -18.5, -16.5);
	fit_q_spectrum.SetParLimits(10, -21.0, -20.0);
	// set limits to sigma
	fit_q_spectrum.SetParLimits(11, 0.1, 10.0);
	fit_q_spectrum.SetNpx(1000);
	hist_q.Fit(&fit_q_spectrum, "R+");

	double q_parameters[12];
	fit_q_spectrum.GetParameters(q_parameters);
	TF1 *fit_q_spectrums[4];
	for (int i = 0; i < 4; ++i) {
		fit_q_spectrums[i] = new TF1(
			TString::Format("fq%d", i), "gaus", -23, -8
		);
		fit_q_spectrums[i]->SetParameters(q_parameters+3*i);
		fit_q_spectrums[i]->SetLineColor(kBlue);
		fit_q_spectrums[i]->SetNpx(200);
		hist_q.GetListOfFunctions()->Add(fit_q_spectrums[i]);
	}

	for (int i = 0; i < 4; ++i) {
		std::cout << q_parameters[i*3] << " " << q_parameters[i*3+1]
			<< " " << q_parameters[i*3+2] << "\n";
	}

	// save
	hist_q.Write();
	c_spec_all.Write();
	c_spec_0.Write();
	c_spec_1.Write();
	c_spec_2.Write();
	c_spec_3.Write();
	// close files
	output_file.Close();
	return 0;
}