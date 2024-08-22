#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <RooRealVar.h>
#include <RooPlot.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooVoigtian.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooCategory.h>
#include <RooSimultaneous.h>

#include "include/defs.h"

#define VOIGT

using namespace ribll;

double inline Resolution0(double x) {
	return -0.649 + 0.290 * log(x);
}


double inline Resolution1(double x) {
	return -0.749 + 0.305 * log(x);
}


double inline Resolution2(double x) {
	return -1.100 + 0.403 * log(x);
}

constexpr double sigma_range = 0.1;

void PlotPdf(
	const RooAbsPdf &pdf,
	RooPlot *frame,
	RooCategory &category,
	const char *label,
	const RooDataHist &hist,
	const RooAbsPdf &component
) {
	pdf.plotOn(
		frame,
		RooFit::Slice(category, label),
		RooFit::ProjWData(category, hist),
		RooFit::LineColor(kBlue),
		RooFit::Components(component)
	);
}

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
	ipt->SetBranchAddress("stateless_excited_energy", stateless_excited_energy);


	// output file name
	TString output_file_name = TString::Format(
		"%s%sthreebody-roofit.root", kGenerateDataPath, kSpectrumDir
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
	TTree *fit_tree[3];
	double fit_ex[3];
	for (int i = 0; i < 3; ++i) {
		fit_tree[i] = new TTree("tree", "tree for RooFit");
		fit_tree[i]->Branch("ex", fit_ex+i, "ex/D");
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
			fit_ex[0] = stateless_excited_energy[0][3];
			fit_tree[0]->Fill();
		}
		if (q[3] < -14.5 && q[3] > -16) {
			hist_ex[1].Fill(stateless_excited_energy[1][3]);
			fit_ex[1] = stateless_excited_energy[1][3];
			fit_tree[1]->Fill();
		}
		if (q[3] < -17.5 && q[3] > -20.0) {
			hist_ex[2].Fill(stateless_excited_energy[2][3]);
			fit_ex[2] = stateless_excited_energy[2][3];
			fit_tree[2]->Fill();
		}
	}

	for (int i = 0; i < 3; ++i) hist_ex[i].Write(TString::Format("ho%d", i));


	// RooFit
	// get data
	RooRealVar x("ex", "excited energy", 12.0, 27.0);
	RooDataHist hist0("dh0", "excited energy", x, hist_ex);
	// RooDataSet set0("s0", "excited energy", fit_tree[0], RooArgSet(x));
	// get gaussian function
	// gaussian at 13.6 MeV
	RooRealVar mean136("mean136", "mean", 13.65, 13.6, 13.7);
#ifdef VOIGT
	// double sigma136_value = Resolution0(13.6);
	// RooRealVar sigma136(
	// 	"sigma136", "sigma",
	// 	sigma136_value,
	// 	sigma136_value-sigma_range,
	// 	sigma136_value+sigma_range
	// );
	// RooRealVar gamma136("gamma136", "gamma", 0.1, 0.0, 0.5);
	// RooVoigtian voigt136("voigt136", "voigt", x, mean136, gamma136, sigma136);
	RooRealVar sigma136("sigma136", "sigma", 0.15, 0.1, 0.2);
	RooGaussian gaus136("gaus136", "gaus", x, mean136, sigma136);
#else
	RooRealVar sigma136("sigma136", "sigma", 0.15, 0.1, 0.2);
	RooGaussian gaus136("gaus136", "gaus", x, mean136, sigma136);
#endif
	// gaussian at 14.9 MeV
	RooRealVar mean149("mean149", "mean", 14.9, 14.7, 15.0);
#ifdef VOIGT
	double sigma149_value = Resolution0(14.9);
	RooRealVar sigma149(
		"sigma149", "sigma",
		sigma149_value,
		sigma149_value-sigma_range,
		sigma149_value+sigma_range
	);
	RooRealVar gamma149("gamma149", "gamma", 0.1, 0.0, 0.5);
	RooVoigtian voigt149("voigt149", "voigt", x, mean149, gamma149, sigma149);
#else
	RooRealVar sigma149("sigma149", "sigma", 0.2, 0.05, 0.5);
	RooGaussian gaus149("gaus149", "gaus", x, mean149, sigma149);
#endif
	// gaussian at 15.4 MeV
	RooRealVar mean154("mean154", "mean", 15.4, 15.3, 15.8);
#ifdef VOIGT
	double sigma154_value = Resolution0(15.4);
	RooRealVar sigma154(
		"sigma154", "sigma",
		sigma154_value,
		sigma154_value-sigma_range,
		sigma154_value+sigma_range
	);
	RooRealVar gamma154("gamma154", "gamma", 0.1, 0.0, 0.5);
	RooVoigtian voigt154("voigt154", "voigt", x, mean154, gamma154, sigma154);
#else
	RooRealVar sigma154("sigma154", "sigma", 0.2, 0.05, 0.5);
	RooGaussian gaus154("gaus154", "gaus", x, mean154, sigma154);
#endif
	// gaussian at 16.5 MeV
	RooRealVar mean165("mean165", "mean", 16.4, 16.2, 16.6);
#ifdef VOIGT
	double sigma165_value = Resolution0(16.5);
	RooRealVar sigma165(
		"sigma165", "sigma",
		sigma165_value,
		sigma165_value-sigma_range,
		sigma165_value+sigma_range
	);
	RooRealVar gamma165("gamma165", "gamma", 0.1, 0.0, 0.5);
	RooVoigtian voigt165("voigt165", "voigt", x, mean165, gamma165, sigma165);
#else
	RooRealVar sigma165("sigma165", "sigma", 0.3, 0.05, 0.5);
	RooGaussian gaus165("gaus165", "gaus", x, mean165, sigma165);
#endif
	// gaussian at 17.2 MeV
	RooRealVar mean172("mean172", "mean", 17.2, 17.0, 17.3);
#ifdef VOIGT
	double sigma172_value = Resolution0(17.2);
	RooRealVar sigma172(
		"sigma172", "sigma",
		sigma172_value,
		sigma172_value-sigma_range,
		sigma172_value+sigma_range
	);
	RooRealVar gamma172("gamma172", "gamma", 0.2, 0.1, 0.5);
	RooVoigtian voigt172("voigt172", "voigt", x, mean172, gamma172, sigma172);
#else
	RooRealVar sigma172("sigma172", "sigma", 0.1, 0.05, 0.5);
	RooGaussian gaus172("gaus172", "gaus", x, mean172, sigma172);
#endif
	// gaussian at 17.7 MeV
	RooRealVar mean177("mean177", "mean", 17.7, 17.4, 17.8);
#ifdef VOIGT
	double sigma177_value = Resolution0(17.7);
	RooRealVar sigma177(
		"sigma177", "sigma",
		sigma177_value,
		sigma177_value-sigma_range,
		sigma177_value+sigma_range
	);
	RooRealVar gamma177("gamma177", "gamma", 0.1, 0.05, 0.5);
	RooVoigtian voigt177("voigt177", "voigt", x, mean177, gamma177, sigma177);
#else
	RooRealVar sigma177("sigma177", "sigma", 0.4, 0.05, 0.5);
	RooGaussian gaus177("gaus177", "gaus", x, mean177, sigma177);
#endif
	// gaussian at 18.5 MeV
	RooRealVar mean185("mean185", "mean", 18.5, 18.0, 18.7);
#ifdef VOIGT
	double sigma185_value = Resolution0(18.5);
	RooRealVar sigma185(
		"sigma185", "sigma",
		sigma185_value,
		sigma185_value-sigma_range,
		sigma185_value+sigma_range
	);
	RooRealVar gamma185("gamma185", "gamma", 0.2, 0.1, 0.5);
	RooVoigtian voigt185("voigt185", "voigt", x, mean185, gamma185, sigma185);
#else
	RooRealVar sigma185("sigma185", "sigma", 0.4, 0.05, 0.5);
	RooGaussian gaus185("gaus185", "gaus", x, mean185, sigma185);
#endif
	// gaussian at 19.3 MeV
	RooRealVar mean193("mean193", "mean", 19.3, 19.2, 19.6);
#ifdef VOIGT
	double sigma193_value = Resolution0(19.3);
	RooRealVar sigma193(
		"sigma193", "sigma",
		sigma193_value,
		sigma193_value-sigma_range,
		sigma193_value+sigma_range
	);
	RooRealVar gamma193("gamma193", "gamma", 0.1, 0.0, 0.5);
	RooVoigtian voigt193("voigt193", "voigt", x, mean193, gamma193, sigma193);
#else
	RooRealVar sigma193("sigma193", "sigma", 0.2, 0.05, 0.5);
	RooGaussian gaus193("gaus193", "gaus", x, mean193, sigma193);
#endif
	// gaussian at 20.1 MeV
	RooRealVar mean201("mean201", "mean", 20.1, 19.8, 20.4);
#ifdef VOIGT
	double sigma201_value = Resolution0(20.1);
	RooRealVar sigma201(
		"sigma201", "sigma",
		sigma201_value,
		sigma201_value-sigma_range,
		sigma201_value+sigma_range
	);
	RooRealVar gamma201("gamma201", "gamma", 0.1, 0.0, 0.5);
	RooVoigtian voigt201("voigt201", "voigt", x, mean201, gamma201, sigma201);
#else
	RooRealVar sigma201("sigma201", "sigma", 0.5, 0.05, 0.5);
	RooGaussian gaus201("gaus201", "gaus", x, mean201, sigma201);
#endif
	// gaussian at 21.4 MeV
	RooRealVar mean214("mean214", "mean", 21.0, 20.5, 21.7);
#ifdef VOIGT
	double sigma214_value = Resolution0(21.0);
	RooRealVar sigma214(
		"sigma214", "sigma",
		sigma214_value,
		sigma214_value-sigma_range,
		sigma214_value+sigma_range
	);
	RooRealVar gamma214("gamma214", "gamma", 0.1, 0.0, 0.5);
	RooVoigtian voigt214("voigt214", "voigt", x, mean214, gamma214, sigma214);
#else
	RooRealVar sigma214("sigma214", "sigma", 0.4, 0.05, 1.0);
	RooGaussian gaus214("gaus214", "gaus", x, mean214, sigma214);
#endif
	// gaussian at 22.4 MeV
	RooRealVar mean224("mean224", "mean", 22.4, 22.0, 23.0);
#ifdef VOIGT
	double sigma224_value = Resolution0(22.4);
	RooRealVar sigma224(
		"sigma224", "sigma",
		sigma224_value,
		sigma224_value-sigma_range,
		sigma224_value+sigma_range
	);
	RooRealVar gamma224("gamma224", "gamma", 0.1, 0.0, 0.5);
	RooVoigtian voigt224("voigt224", "voigt", x, mean224, gamma224, sigma224);
#else
	RooRealVar sigma224("sigma224", "sigma", 0.4, 0.05, 0.6);
	RooGaussian gaus224("gaus224", "gaus", x, mean224, sigma224);
#endif
	// gaussian at 23.8 MeV
	RooRealVar mean238("mean238", "mean", 23.8, 23.0, 24.0);
#ifdef VOIGT
	double sigma238_value = Resolution0(23.8);
	RooRealVar sigma238(
		"sigma238", "sigma",
		sigma238_value,
		sigma238_value-sigma_range,
		sigma238_value+sigma_range
	);
	RooRealVar gamma238("gamma238", "gamma", 0.1, 0.0, 0.5);
	RooVoigtian voigt238("voigt238", "voigt", x, mean238, gamma238, sigma238);
#else
	RooRealVar sigma238("sigma238", "sigma", 0.4, 0.05, 1.0);
	RooGaussian gaus238("gaus238", "gaus", x, mean238, sigma238);
#endif
	// gaussian at 24.5 MeV
	RooRealVar mean245("mean245", "mean", 24.5, 24.2, 24.7);
#ifdef VOIGT
	double sigma245_value = Resolution0(24.5);
	RooRealVar sigma245(
		"sigma245", "sigma",
		sigma245_value,
		sigma245_value-sigma_range,
		sigma245_value+sigma_range
	);
	RooRealVar gamma245("gamma245", "gamma", 0.1, 0.0, 0.5);
	RooVoigtian voigt245("voigt245", "voigt", x, mean245, gamma245, sigma245);
#else
	RooRealVar sigma245("sigma245", "sigma", 0.2, 0.05, 0.6);
	RooGaussian gaus245("gaus245", "gaus", x, mean245, sigma245);
#endif
	// gaussian at 25.5 MeV
	RooRealVar mean259("mean259", "mean", 25.9, 25, 26.4);
#ifdef VOIGT
	double sigma259_value = Resolution0(25.9);
	RooRealVar sigma259(
		"sigma259", "sigma",
		sigma259_value,
		sigma259_value-sigma_range,
		sigma259_value+sigma_range
	);
	RooRealVar gamma259("gamma259", "gamma", 0.1, 0.0, 0.5);
	RooVoigtian voigt259("voigt259", "voigt", x, mean259, gamma259, sigma259);
#else
	RooRealVar sigma259("sigma259", "sigma", 0.7, 0.05, 1.0);
	RooGaussian gaus259("gaus259", "gaus", x, mean259, sigma259);
#endif
	// summary
	RooAbsPdf* fit_ex0[16] = {
#ifdef VOIGT
		&gaus136, &voigt149, &voigt154, &voigt165,
		&voigt172, &voigt177, &voigt185, &voigt193,
		&voigt201, &voigt214, &voigt224, &voigt238,
		&voigt245, &voigt259, nullptr, nullptr
#else
		&gaus136, &gaus149, &gaus154, &gaus165,
		&gaus172, &gaus177, &gaus185, &gaus193,
		&gaus201, &gaus214, &gaus224, &gaus238,
		&gaus245, &gaus259, nullptr, nullptr
#endif
	};
	RooRealVar* frac0[16];
	RooAbsPdf* model0[16];
	// add gaussian
	// number of peaks
	size_t n0 = 13;
	model0[0] = fit_ex0[0];
	for (size_t i = 1; i < n0; ++i) {
		frac0[i-1] = new RooRealVar(
			TString::Format("frac0_%ld", i), "frac",
			1.0-1.0/i, 0.0, 1.0
		);
		model0[i] = new RooAddPdf(
			TString::Format("model0_%ld", i), "model",
			RooArgList(*model0[i-1], *fit_ex0[i]),
			*frac0[i-1]
		);
	}


	// 10Be 3.5 MeV
	// RooRealVar x("ex", "ex1", 12.0, 27.0);
	RooDataHist hist1("dh1", "excited energy", x, hist_ex+1);
	// RooDataSet set1("s1", "excited energy", fit_tree[1], RooArgSet(x1));
	// gaussian at 17.2 MeV
#ifdef VOIGT
	double sigma172_1_value = Resolution1(17.2);
	RooRealVar sigma172_1(
		"sigma172_1", "sigma",
		sigma172_1_value,
		sigma172_1_value-sigma_range,
		sigma172_1_value+sigma_range
	);
	RooVoigtian voigt172_1(
		"voigt172_1", "voigt", x, mean172, gamma172, sigma172_1
	);
	// RooRealVar sigma172_1("sigma172_1", "sigma", 0.15, 0.1, 0.2);
	// RooGaussian gaus172_1("gaus172_1", "gaus", x, mean172, sigma172_1);
#else
	RooRealVar sigma172_1("sigma172_1", "sigma", 0.15, 0.1, 0.2);
	RooGaussian gaus172_1("gaus172_1", "gaus", x, mean172, sigma172_1);
#endif
	// gaussian at 17.7-17.9 MeV
#ifdef VOIGT
	double sigma177_1_value = Resolution1(17.7);
	RooRealVar sigma177_1(
		"sigma177_1", "sigma",
		sigma177_1_value,
		sigma177_1_value-sigma_range,
		sigma177_1_value+sigma_range
	);
	RooVoigtian voigt177_1(
		"voigt177_1", "voigt", x, mean177, gamma177, sigma177_1
	);
#else
	RooRealVar sigma177_1("sigma177_1", "sigma", 0.2, 0.05, 1.0);
	RooGaussian gaus177_1("gaus177_1", "gaus", x, mean177, sigma177_1);
#endif
	// gaussian at 18.5 MeV
#ifdef VOIGT
	double sigma185_1_value = Resolution1(18.5);
	RooRealVar sigma185_1(
		"sigma185_1", "sigma",
		sigma185_1_value,
		sigma185_1_value-sigma_range,
		sigma185_1_value+sigma_range
	);
	RooVoigtian voigt185_1(
		"voigt185_1", "voigt", x, mean185, gamma185, sigma185_1
	);
#else
	RooRealVar sigma185_1("sigma185_1", "sigma", 0.2, 0.15, 0.5);
	RooGaussian gaus185_1("gaus185_1", "gaus", x, mean185, sigma185_1);
#endif
	// gaussian at 19.3 MeV
#ifdef VOIGT
	double sigma193_1_value = Resolution1(19.3);
	RooRealVar sigma193_1(
		"sigma193_1", "sigma",
		sigma193_1_value,
		sigma193_1_value-sigma_range,
		sigma193_1_value+sigma_range
	);
	RooVoigtian voigt193_1(
		"voigt193_1", "voigt", x, mean193, gamma193, sigma193_1
	);
#else
	RooRealVar sigma193_1("sigma193_1", "sigma", 0.2, 0.05, 0.5);
	RooGaussian gaus193_1("gaus193_1", "gaus", x, mean193, sigma193_1);
#endif
	// gaussian at 20.1 MeV
#ifdef VOIGT
	double sigma201_1_value = Resolution1(20.1);
	RooRealVar sigma201_1(
		"sigma201_1", "sigma",
		sigma201_1_value,
		sigma201_1_value-sigma_range,
		sigma201_1_value+sigma_range
	);
	RooVoigtian voigt201_1(
		"voigt201_1", "voigt", x, mean201, gamma201, sigma201_1
	);
#else
	RooRealVar sigma201_1("sigma201_1", "sigma", 0.2, 0.05, 0.5);
	RooGaussian gaus201_1("gaus201_1", "gaus", x, mean201, sigma201_1);
#endif
	// gaussian at 21.4 MeV
#ifdef VOIGT
	double sigma214_1_value = Resolution1(21.4);
	RooRealVar sigma214_1(
		"sigma214_1", "sigma",
		sigma214_1_value,
		sigma214_1_value-sigma_range,
		sigma214_1_value+sigma_range
	);
	RooVoigtian voigt214_1(
		"voigt214_1", "voigt", x, mean214, gamma214, sigma214_1
	);
#else
	RooRealVar sigma214_1("sigma214_1", "sigma", 0.3, 0.05, 1.0);
	RooGaussian gaus214_1("gaus214_1", "gaus", x, mean214, sigma214_1);
#endif
	// gaussian at 22.4 MeV
#ifdef VOIGT
	double sigma224_1_value = Resolution1(22.4);
	RooRealVar sigma224_1(
		"sigma224_1", "sigma",
		sigma224_1_value,
		sigma224_1_value-sigma_range,
		sigma224_1_value+sigma_range
	);
	RooVoigtian voigt224_1(
		"voigt224_1", "voigt", x, mean224, gamma224, sigma224_1
	);
#else
	RooRealVar sigma224_1("sigma224_1", "sigma", 0.2, 0.05, 0.6);
	RooGaussian gaus224_1("gaus224_1", "gaus", x, mean224, sigma224_1);
#endif
	// gaussian at 23.5 MeV
	RooRealVar mean235("mean235", "mean", 23.5, 23.2, 23.9);
#ifdef VOIGT
	double sigma235_value = Resolution1(23.5);
	RooRealVar sigma235(
		"sigma235", "sigma",
		sigma235_value,
		sigma235_value-sigma_range,
		sigma235_value+sigma_range
	);
	RooRealVar gamma235("gamma235", "gamma", 0.1, 0.0, 0.5);
	RooVoigtian voigt235(
		"voigt235", "voigt", x, mean235, gamma235, sigma235
	);
#else
	RooRealVar sigma235("sigma235", "sigma", 0.2, 0.05, 1.0);
	RooGaussian gaus235("gaus235", "gaus", x, mean235, sigma235);
#endif
	// gaussian at 24.5 MeV
#ifdef VOIGT
	double sigma245_1_value = Resolution1(24.5);
	RooRealVar sigma245_1(
		"sigma245_1", "sigma",
		sigma245_1_value,
		sigma245_1_value-sigma_range,
		sigma245_1_value+sigma_range
	);
	RooVoigtian voigt245_1(
		"voigt245_1", "voigt", x, mean245, gamma245, sigma245_1
	);
#else
	RooRealVar sigma245_1("sigma245_1", "sigma", 0.3, 0.05, 1.0);
	RooGaussian gaus245_1("gaus245_1", "gaus", x, mean245, sigma245_1);
#endif
	// gaussian at 25.9 MeV
#ifdef VOIGT
	double sigma259_1_value = Resolution1(25.9);
	RooRealVar sigma259_1(
		"sigma259_1", "sigma",
		sigma259_1_value,
		sigma259_1_value-sigma_range,
		sigma259_1_value+sigma_range
	);
	RooVoigtian voigt259_1(
		"voigt259_1", "voigt", x, mean259, gamma259, sigma259_1
	);
#else
	RooRealVar sigma259_1("sigma259_1", "sigma", 0.2, 0.05, 1.0);
	RooGaussian gaus259_1("gaus259_1", "gaus", x, mean259, sigma259_1);
#endif
	// summary
	RooAbsPdf* fit_ex1[16] = {
#ifdef VOIGT
		&voigt172_1, &voigt177_1, &voigt185_1, &voigt193_1,
		&voigt201_1, &voigt214_1, &voigt224_1, &voigt235,
		&voigt245_1, &voigt259_1, nullptr, nullptr,
		nullptr, nullptr, nullptr, nullptr
#else
		&gaus172_1, &gaus177_1, &gaus185_1, &gaus193_1,
		&gaus201_1, &gaus214_1, &gaus224_1, &gaus235,
		&gaus245_1, &gaus259_1, nullptr, nullptr,
		nullptr, nullptr, nullptr, nullptr
#endif
	};
	RooRealVar* frac1[16];
	RooAbsPdf* model1[16];
	// add gaussian
	// number of peaks
	size_t n1 = 9;
	model1[0] = fit_ex1[0];
	for (size_t i = 1; i < n1; ++i) {
		frac1[i-1] = new RooRealVar(
			TString::Format("frac1_%ld", i), "frac",
			1.0-1.0/i, 0.0, 1.0
		);
		model1[i] = new RooAddPdf(
			TString::Format("model1_%ld", i), "model",
			RooArgList(*model1[i-1], *fit_ex1[i]),
			*frac1[i-1]
		);
	}

	// 10Be 6 MeV
	RooDataHist hist2("dh2", "excited energy", RooArgList(x), hist_ex+2);
	RooDataSet set2("s2", "excited energy", fit_tree[2], RooArgSet(x));
	// gaussian at 20.1 MeV
#ifdef VOIGT
	double sigma201_2_value = Resolution2(20.1);
	RooRealVar sigma201_2(
		"sigma201_2", "sigma",
		sigma201_2_value,
		sigma201_2_value-sigma_range,
		sigma201_2_value+sigma_range
	);
	RooVoigtian voigt201_2(
		"voigt201_2", "voigt", x, mean201, gamma201, sigma201_2
	);
#else
	RooRealVar sigma201_2("sigma201_2", "sigma", 0.4, 0.05, 1.0);
	RooGaussian gaus201_2("gaus201_2", "gaus", x, mean201, sigma201_2);
#endif
	// gaussian at 21.4 MeV
#ifdef VOIGT
	double sigma214_2_value = Resolution2(21.4);
	RooRealVar sigma214_2(
		"sigma214_2", "sigma",
		sigma214_2_value,
		sigma214_2_value-sigma_range,
		sigma214_2_value+sigma_range
	);
	RooVoigtian voigt214_2(
		"voigt214_2", "voigt", x, mean214, gamma214, sigma214_2
	);
#else
	RooRealVar sigma214_2("sigma214_2", "sigma", 0.5, 0.05, 1.0);
	RooGaussian gaus214_2("gaus214_2", "gaus", x, mean214, sigma214_2);
#endif
	// gaussian at 22.4 MeV
#ifdef VOIGT
	double sigma224_2_value = Resolution2(22.4);
	RooRealVar sigma224_2(
		"sigma224_2", "sigma",
		sigma224_2_value,
		sigma224_2_value-sigma_range,
		sigma224_2_value+sigma_range
	);
	RooVoigtian voigt224_2(
		"voigt224_2", "voigt", x, mean224, gamma224, sigma224_2
	);
#else
	RooRealVar sigma224_2("sigma224_2", "sigma", 0.5, 0.05, 1.0);
	RooGaussian gaus224_2("gaus224_2", "gaus", x, mean224, sigma224_2);
#endif
	// gaussian at 23.5 MeV
#ifdef VOIGT
	double sigma235_2_value = Resolution2(23.5);
	RooRealVar sigma235_2(
		"sigma235_2", "sigma",
		sigma235_2_value,
		sigma235_2_value-sigma_range,
		sigma235_2_value+sigma_range
	);
	RooVoigtian voigt235_2(
		"voigt235_2", "voigt", x, mean235, gamma235, sigma235_2
	);
#else
	RooRealVar sigma235_2("sigma235_2", "sigma", 0.7, 0.05, 1.0);
	RooGaussian gaus235_2("gaus235_2", "gaus", x, mean235, sigma235_2);
#endif
	// gaussian at 24.5 MeV
#ifdef VOIGT
	double sigma245_2_value = Resolution2(24.5);
	RooRealVar sigma245_2(
		"sigma245_2", "sigma",
		sigma245_2_value,
		sigma245_2_value-sigma_range,
		sigma245_2_value+sigma_range
	);
	RooVoigtian voigt245_2(
		"voigt245_2", "voigt", x, mean245, gamma245, sigma245_2
	);
#else
	RooRealVar sigma245_2("sigma245_2", "sigma", 0.2, 0.05, 1.0);
	RooGaussian gaus245_2("gaus245_2", "gaus", x, mean245, sigma245_2);
#endif
	// gaussian at 25.9 MeV
#ifdef VOIGT
	double sigma259_2_value = Resolution2(25.9);
	RooRealVar sigma259_2(
		"sigma259_2", "sigma",
		sigma259_2_value,
		sigma259_2_value-sigma_range,
		sigma259_2_value+sigma_range
	);
	RooVoigtian voigt259_2(
		"voigt259_2", "voigt", x, mean259, gamma259, sigma259_2
	);
#else
	RooRealVar sigma259_2("sigma259_2", "sigma", 0.6, 0.05, 1.0);
	RooGaussian gaus259_2("gaus259_2", "gaus", x, mean259, sigma259_2);
#endif
	// summary
	RooAbsPdf* fit_ex2[16] = {
#ifdef VOIGT
		&voigt201_2, &voigt214_2, &voigt224_2, &voigt235_2,
		&voigt245_2, &voigt259_2, nullptr, nullptr,
		nullptr, nullptr, nullptr, nullptr,
		nullptr, nullptr, nullptr, nullptr
#else
		&gaus201_2, &gaus214_2, &gaus224_2, &gaus235_2,
		&gaus245_2, &gaus259_2, nullptr, nullptr,
		nullptr, nullptr, nullptr, nullptr,
		nullptr, nullptr, nullptr, nullptr
#endif
	};
	RooRealVar* frac2[16];
	RooAbsPdf* model2[16];
	// add gaussian
	// number of peaks
	size_t n2 = 6;
	model2[0] = fit_ex2[0];
	for (size_t i = 1; i < n2; ++i) {
		frac2[i-1] = new RooRealVar(
			TString::Format("frac2_%ld", i), "frac",
			1.0-1.0/i, 0.0, 1.0
		);
		model2[i] = new RooAddPdf(
			TString::Format("model2_%ld", i), "model",
			RooArgList(*model2[i-1], *fit_ex2[i]),
			*frac2[i-1]
		);
	}


	// simultaneous fitting
	// category
	RooCategory state("state", "state");
	state.defineType("0");
	state.defineType("1");
	state.defineType("2");
	// combined data
	RooDataHist combined_hist(
		"combined", "combined data",
		x,
		RooFit::Index(state),
		RooFit::Import("0", hist0),
		RooFit::Import("1", hist1),
		RooFit::Import("2", hist2)
	);
	combined_hist.Print();
	// simultaneous P.D.F.
	RooSimultaneous simultaneous_model(
		"sim_model", "simultaneous pdf",
		{{"0", model0[n0-1]}, {"1", model1[n1-1]}, {"2", model2[n2-1]}},
		state
	);
	simultaneous_model.Print();
	std::unique_ptr<RooFitResult> fit_result{
		simultaneous_model.fitTo(
			combined_hist,
			RooFit::PrintLevel(-1),
			RooFit::Save()
		)
	};
	fit_result->Print();

	// plot
	// frame 0
	RooPlot *frame0 = x.frame();
	combined_hist.plotOn(
		frame0,
		RooFit::Cut("state==state::0"),
		RooFit::DrawOption("B"),
		RooFit::FillColor(kCyan),
		RooFit::DataError(RooAbsData::None)
	);
	simultaneous_model.plotOn(
		frame0,
		RooFit::Slice(state, "0"),
		RooFit::ProjWData(state, combined_hist),
		RooFit::LineColor(kRed)
	);
	for (size_t i = 0; i < n0; ++i) {
		PlotPdf(
			simultaneous_model, frame0,
			state, "0",
			combined_hist, *fit_ex0[i]
		);
	}

	// frame 1
	RooPlot *frame1 = x.frame();
	combined_hist.plotOn(
		frame1,
		RooFit::Cut("state==state::1"),
		RooFit::DrawOption("B"),
		RooFit::FillColor(kCyan),
		RooFit::DataError(RooAbsData::None)
	);
	simultaneous_model.plotOn(
		frame1,
		RooFit::Slice(state, "1"),
		RooFit::ProjWData(state, combined_hist),
		RooFit::LineColor(kRed)
	);
	for (size_t i = 0; i < n1; ++i) {
		PlotPdf(
			simultaneous_model, frame1,
			state, "1",
			combined_hist, *fit_ex1[i]
		);
	}

	// frame 2
	RooPlot *frame2 = x.frame();
	combined_hist.plotOn(
		frame2,
		RooFit::Cut("state==state::2"),
		RooFit::DrawOption("B"),
		RooFit::FillColor(kCyan),
		RooFit::DataError(RooAbsData::None)
	);
	simultaneous_model.plotOn(
		frame2,
		RooFit::Slice(state, "2"),
		RooFit::ProjWData(state, combined_hist),
		RooFit::LineColor(kRed)
	);
	for (size_t i = 0; i < n2; ++i) {
		PlotPdf(
			simultaneous_model, frame2,
			state, "2",
			combined_hist, *fit_ex2[i]
		);
	}

	// close files
	for (int i = 0; i < 3; ++i) hist_ex[i].Write();
	frame0->Write("rh0");
	frame1->Write("rh1");
	frame2->Write("rh2");

	// save
	gStyle->SetPadBorderMode(0);
	TCanvas *c1 = new TCanvas("c1", "", 1920, 1080);
	// build pads
	TPad *pads[3];
	for (int i = 0; i < 3; ++i) {
		pads[i] = new TPad(
			TString::Format("pad%d", i), "",
			0, (7-i*3)/11.0, 1, (10-i*3)/11.0
		);
		pads[i]->Draw();
	}
	pads[0]->SetBottomMargin(0);
	pads[1]->SetTopMargin(0);
	pads[1]->SetBottomMargin(0);
	pads[2]->SetTopMargin(0);
	pads[0]->cd();
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	frame0->GetYaxis()->SetLabelSize(0.15);
	frame0->Draw();
	pads[1]->cd();
	frame1->Draw();
	frame1->GetYaxis()->SetLabelSize(0.15);
	pads[2]->cd();
	frame2->Draw();
	frame2->GetXaxis()->SetLabelSize(0.12);
	frame2->GetYaxis()->SetLabelSize(0.15);
	c1->SaveAs(
		TString::Format("%simage/rootfit-result.png", kGenerateDataPath)
	);
	c1->Write("c1");

	opf.Close();
	ipf.Close();

	// clean
	for (size_t i = 1; i < n0; ++i) {
		delete frac0[i-1];
		delete model0[i];
	}
	for (size_t i = 1; i < n1; ++i) {
		delete frac1[i-1];
		delete model1[i];
	}
	for (size_t i = 1; i < n2; ++i) {
		delete frac2[i-1];
		delete model2[i];
	}
	return 0;
}