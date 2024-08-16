#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>

#include "include/defs.h"

using namespace ribll;

class GausFunctor {
public:
	GausFunctor(int n): n_(n) {}

	double operator()(double *x, double *par) const {
		double result = 0.0;
		for (int i = 0; i < n_; ++i) {
			result += par[i*3]*exp(-0.5*pow((x[0]-par[i*3+1])/par[i*3+2], 2.0));
		}
		return result;
	}

private:
	int n_;
};

int main() {
	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody-10Be-dv2-2.root", kGenerateDataPath, kSpectrumDir
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// input data
	int ppac_flag, taf_flag, target_flag, bind, hole;
	double q[4];
	int be_state[4];
	double excited_energy_target[4];
	double stateless_excited_energy[3][4];
	// setup input branches
	ipt->SetBranchAddress("ppac_flag", &ppac_flag);
	ipt->SetBranchAddress("taf_flag", &taf_flag);
	ipt->SetBranchAddress("target_flag", &target_flag);
	ipt->SetBranchAddress("bind", &bind);
	ipt->SetBranchAddress("hole", &hole);
	ipt->SetBranchAddress("q", q);
	ipt->SetBranchAddress("be_state", be_state);
	ipt->SetBranchAddress("excited_energy_target", excited_energy_target);
	ipt->SetBranchAddress("stateless_excited_energy", stateless_excited_energy);;


	// output file name
	TString output_file_name = TString::Format(
		"%s%sthreebody-fit.root", kGenerateDataPath, kSpectrumDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	TH1F hist_ex[3];
	for (int i = 0; i < 3; ++i) {
		hist_ex[i] = TH1F(
			TString::Format("h%d", i), "excited energy", 50, 12, 27
		);
		hist_ex[i].SetLineColor(kBlack);
	}

	// loop
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);
		if ((ppac_flag & 1) != 1) continue;
		if (taf_flag != 0) continue;
		if ((target_flag & 1) == 0) continue;
		if (bind != 0) continue;
		if (hole != 0) continue;
		if (q[3] < -10 && q[3] > -14) {
			hist_ex[0].Fill(stateless_excited_energy[0][3]);
		}
		if (q[3] < -14.5 && q[3] > -16) {
			hist_ex[1].Fill(stateless_excited_energy[1][3]);
		}
		if (q[3] < -17.5 && q[3] > -20.0) {
			hist_ex[2].Fill(stateless_excited_energy[2][3]);
		}
	}

	for (int i = 0; i < 3; ++i) hist_ex[i].Write(TString::Format("ho%d", i));


	// fit ground state
	const int n0 = 8;
	GausFunctor fit_s0(n0);
	double init_param0[n0*3] = {
		25.0, 15.5, 0.5,
		20.0, 16.4, 0.5,
		20.0, 18.0, 0.5,
		20.0, 19.5, 0.5,
		20.0, 21.0, 0.5,
		15.0, 22.0, 0.5,
		10.0, 23.5, 0.2,
		10.0, 24.5, 0.2
	};
	TF1 *fs0 = new TF1("fs0", fit_s0, 14.0, 25.0, n0*3);
	fs0->SetNpx(1000);
	fs0->SetParameters(init_param0);
	fs0->SetParLimits(1, 15.3, 15.7);
	fs0->SetParLimits(2, 0.02, 1.0);
	fs0->SetParLimits(4, 16.2, 16.6);
	fs0->SetParLimits(5, 0.05, 0.7);
	fs0->SetParLimits(7, 17.0, 18.3);
	fs0->SetParLimits(8, 0.05, 1.0);
	fs0->SetParLimits(10, 19.0, 20.0);
	fs0->SetParLimits(11, 0.05, 1.0);
	fs0->SetParLimits(13, 20.5, 21.5);
	fs0->SetParLimits(14, 0.05, 1.0);
	fs0->SetParLimits(16, 22.0, 23.0);
	fs0->SetParLimits(17, 0.05, 1.0);
	fs0->SetParLimits(19, 23.0, 24.0);
	fs0->SetParLimits(20, 0.05, 1.0);
	fs0->SetParLimits(22, 24.0, 25.0);
	fs0->SetParLimits(23, 0.05, 1.0);
	hist_ex[0].Fit(fs0, "R+");
	for (int i = 0; i < n0; ++i) {
		TF1 *f1 = new TF1(TString::Format("fs0p%d", i), "gaus", 14, 27);
		f1->SetParameters(fs0->GetParameters()+i*3);
		f1->SetLineColor(kBlue);
		f1->SetNpx(1000);
		hist_ex[0].GetListOfFunctions()->Add(f1);
	}


	// fit 3.5MeV state
	const int n1 = 9;
	GausFunctor fit_s1(n1);
	double init_param1[n1*3] = {
		10.0, 17.3, 0.5,
		20.0, 18.5, 0.5,
		40.0, 19.5, 0.5,
		40.0, 20.2, 0.5,
		40.0, 21.2, 0.5,
		40.0, 22.2, 0.5,
		40.0, 23.5, 0.5,
		10.0, 24.5, 0.5,
		20.0, 26.0, 0.5
	};
	TF1 *fs1 = new TF1("fs1", fit_s1, 16.0, 26.5, n1*3);
	fs1->SetNpx(1000);
	fs1->SetParameters(init_param1);
	fs1->SetParLimits(1, 17.0, 17.8);
	fs1->SetParLimits(2, 0.02, 1.0);
	fs1->SetParLimits(4, 18.2, 18.6);
	fs1->SetParLimits(5, 0.05, 1.0);
	fs1->SetParLimits(7, 19.0, 19.6);
	fs1->SetParLimits(8, 0.05, 1.0);
	fs1->SetParLimits(10, 19.7, 20.4);
	fs1->SetParLimits(11, 0.05, 1.0);
	fs1->SetParLimits(13, 20.8, 21.4);
	fs1->SetParLimits(14, 0.05, 1.0);
	fs1->SetParLimits(16, 22.0, 22.5);
	fs1->SetParLimits(17, 0.05, 1.0);
	fs1->SetParLimits(19, 23.0, 23.6);
	fs1->SetParLimits(20, 0.05, 1.0);
	fs1->SetParLimits(22, 24.2, 24.8);
	fs1->SetParLimits(23, 0.05, 1.0);
	fs1->SetParLimits(25, 25.4, 26.4);
	fs1->SetParLimits(26, 0.05, 1.0);
	hist_ex[1].Fit(fs1, "R+");
	for (int i = 0; i < n1; ++i) {
		TF1 *f1 = new TF1(TString::Format("fs1p%d", i), "gaus", 16, 27);
		f1->SetParameters(fs1->GetParameters()+i*3);
		f1->SetLineColor(kBlue);
		f1->SetNpx(1000);
		hist_ex[1].GetListOfFunctions()->Add(f1);
	}


	// fit 6MeV state
	const int n2 = 6;
	GausFunctor fit_s2(n2);
	double init_param2[n2*3] = {
		10.0, 20.3, 0.5,
		20.0, 21.3, 0.5,
		30.0, 22.4, 0.5,
		20.0, 23.5, 0.5,
		20.0, 24.2, 0.5,
		10.0, 26.2, 0.5
	};
	TF1 *fs2 = new TF1("fs2", fit_s2, 18.0, 27.0, n2*3);
	fs2->SetNpx(1000);
	fs2->SetParameters(init_param2);
	fs2->SetParLimits(1, 19.7, 20.5);
	fs2->SetParLimits(2, 0.02, 0.5);
	fs2->SetParLimits(4, 20.6, 21.7);
	fs2->SetParLimits(5, 0.05, 0.7);
	fs2->SetParLimits(7, 21.9, 22.5);
	fs2->SetParLimits(8, 0.05, 1.0);
	fs2->SetParLimits(10, 23.0, 24.0);
	fs2->SetParLimits(11, 0.05, 1.0);
	fs2->SetParLimits(13, 24.0, 25.0);
	fs2->SetParLimits(14, 0.05, 1.0);
	fs2->SetParLimits(16, 25.6, 26.4);
	fs2->SetParLimits(17, 0.05, 1.0);
	hist_ex[2].Fit(fs2, "R+");
	for (int i = 0; i < n2; ++i) {
		TF1 *f1 = new TF1(TString::Format("fs2p%d", i), "gaus", 18, 27);
		f1->SetParameters(fs2->GetParameters()+i*3);
		f1->SetLineColor(kBlue);
		f1->SetNpx(1000);
		hist_ex[2].GetListOfFunctions()->Add(f1);
	}

	// close files
	for (int i = 0; i < 3; ++i) hist_ex[i].Write();
	opf.Close();
	ipf.Close();
	return 0;
}