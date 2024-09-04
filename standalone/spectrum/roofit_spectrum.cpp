#include <iostream>
#include <iomanip>
#include <map>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <RooRealVar.h>
#include <RooRealConstant.h>
#include <RooPlot.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooCategory.h>
#include <RooSimultaneous.h>

#include "include/defs.h"
#include "include/spectrum/asymmetric_voigt.h"


using namespace ribll;


struct VoigtInfo {
	RooRealVar *mean;
	RooRealVar *gamma;
	RooConstVar *sigma;
	double threshold;
};

constexpr double threshold[3] = {12.01, 15.38, 18.19};
constexpr double res_correct = 1.0;
constexpr double sigma_range = 0.0;

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

	// efficiency file name
	TString efficiency_file_name = TString::Format(
		"%s%sefficiency-good-0002.root", kGenerateDataPath, kSimulateDir
	);
	// efficiency file
	TFile efficiency_file(efficiency_file_name, "read");
	// efficiency graphs
	TGraph *g_efficiency[3];
	g_efficiency[0] = (TGraph*)efficiency_file.Get("g0");
	g_efficiency[1] = (TGraph*)efficiency_file.Get("g1");
	g_efficiency[2] = (TGraph*)efficiency_file.Get("g2");

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
	TGraph graph_resolution[3];
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
		if (hole > 0) continue;
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


	// test asymmetric voigt
	RooRealVar test_x("tx", "test x", 12.0, 27.0);
	RooRealVar test_mean("tmean", "test  mean", 16.0);
	RooRealVar test_g("tg", "test g", 1.0);
	RooRealVar test_sigma("tsigma", "test sigma", 0.3);
	AsymmetricVoigtian test_avoigt(
		"tav", "test asym voigt",
		test_x, test_mean, test_g, test_sigma,
		12.0125, g_efficiency[0]
	);
	RooPlot *test_frame = test_x.frame();
	test_avoigt.plotOn(
		test_frame,
		RooFit::LineColor(kBlue)
	);

	// fitting function information
	std::map<std::string, VoigtInfo> info_map[3];

	// RooFit
	// get data
	RooRealVar x("ex", "excited energy", 12.0, 27.0);
	RooDataHist hist0("dh0", "excited energy", x, hist_ex);
	// RooDataSet set0("s0", "excited energy", fit_tree[0], RooArgSet(x));
	// number of peaks
	size_t n0 = 13;
	// mean value
	const double mean0_value[16][3] = {
		{13.65, 13.6, 13.7},	// 13.6
		{14.9, 14.5, 15.0},		// 14.9
		{15.4, 15.3, 15.8},		// 15.4
		{16.4, 16.2, 16.6},		// 16.5
		// {17.2, 17.0, 17.3},		// 17.2
		{17.7, 17.4, 17.8},		// 17.7
		{18.5, 18.0, 18.7},		// 18.5
		{19.3, 19.2, 19.6},		// 19.3
		{20.1, 19.8, 20.4},		// 20.1
		{21.0, 20.5, 21.7},		// 21.4
		{22.4, 22.0, 22.8},		// 22.4
		// {23.0, 22.8, 23.5},		// 23.0
		{23.8, 23.5, 24.0},		// 23.8
		{24.5, 24.2, 24.7},		// 24.5
		{25.9, 25.0, 26.4},		// 25.9
	};
	// g factor value
	const double gamma0_value[16][3] = {
		{0.05, 0.0, 1.0},		// 13.6
		{0.06, 0.0, 1.0},		// 14.9
		{0.1, 0.0, 1.0},		// 15.4
		{0.1, 0.0, 1.0},		// 16.5
		// {0.05, 0.0, 1.0},		// 17.2
		{0.05, 0.0, 1.0},		// 17.7
		{0.05, 0.0, 1.0},		// 18.5
		{0.05, 0.0, 1.0},		// 19.3
		{0.02, 0.0, 1.0},		// 20.1
		{0.05, 0.0, 1.0},		// 21.4
		{0.02, 0.0, 1.0},		// 22.4
		// {0.02, 0.0, 1.0},		// 23.0
		{0.02, 0.0, 1.0},		// 23.8
		{0.02, 0.0, 1.0},		// 24.5
		{0.02, 0.0, 1.0},		// 25.9
	};
	// sigma value
	const double sigma0_value[16] = {
		0.2, 0.24, 0.26, 0.31,	// 13.6, 14.9, 15.4, 16.5
		0.32, 0.36, 0.38,	// 17.2, 17.7, 18.5, 19.3
		0.39, 0.4, 0.41,	// 20.1, 21.4, 22.4, 23.0
		0.42, 0.43, 0.44, 		// 23.8, 24.5, 25.9
	};
	// variable name
	const std::string name0[16] = {
		"136", "149", "154", "165",	// 13.6, 14.9, 15.4, 16.5
		"177", "185", "193",	// 17.2, 17.7, 18.5, 19.3
		"201", "214", "224", // 20.1, 21.4, 22.4, 23.0
		"238", "245", "259",		// 23.8, 24.5, 25.9
	};
	// threshold
	const double threshold0[16] = {
		// 13.6, 14.9, 15.4, 16.5
		threshold[0], threshold[0], threshold[0], threshold[0],
		// 17.2, 17.7, 18.5, 19.3
		threshold[1], threshold[1], threshold[1],
		// 20.1, 21.4, 22.4, 23.0
		threshold[2], threshold[2], threshold[2],
		// 23.8, 24.5, 25.9
		threshold[2], threshold[2], threshold[2]
	};
	// mean variable
	RooRealVar *mean0[16];
	// g factor variable
	RooRealVar *gamma0[16];
	// sigma variable
	RooConstVar *sigma0[16];
	// P.D.F
	RooAbsPdf *pdf0[16];
	// create RooAbsPdf
	for (size_t i = 0; i < n0; ++i) {
		mean0[i] = new RooRealVar(
			("mean"+name0[i]).c_str(), "mean",
			mean0_value[i][0], mean0_value[i][1], mean0_value[i][2]
		);
		gamma0[i] = new RooRealVar(
			("gamma"+name0[i]).c_str(), "gamma",
			gamma0_value[i][0], gamma0_value[i][1], gamma0_value[i][2]
		);
		sigma0[i] = new RooConstVar(
			("sigma"+name0[i]).c_str(), "sigma",
			sigma0_value[i]
		);
		pdf0[i] = new AsymmetricVoigtian(
			("voigt"+name0[i]).c_str(), "voigt",
			x, *mean0[i], *gamma0[i], *sigma0[i],
			threshold0[i], g_efficiency[0]
		);
		VoigtInfo info;
		info.mean = mean0[i];
		info.gamma = gamma0[i];
		info.sigma = sigma0[i];
		info.threshold = threshold0[i];
		info_map[0].insert(std::make_pair(name0[i], info));
	}
	// add fraction
	RooRealVar* frac0[16];
	// add result
	RooAbsPdf* model0[16];
	// add gaussian
	model0[0] = pdf0[0];
	for (size_t i = 1; i < n0; ++i) {
		frac0[i-1] = new RooRealVar(
			TString::Format("frac0_%ld", i), "frac",
			1.0-1.0/i, 0.0, 1.0
		);
		model0[i] = new RooAddPdf(
			TString::Format("model0_%ld", i), "model",
			RooArgList(*model0[i-1], *pdf0[i]),
			*frac0[i-1]
		);
	}


	// 10Be 3.5 MeV
	// RooRealVar x("ex", "ex1", 12.0, 27.0);
	RooDataHist hist1("dh1", "excited energy", x, hist_ex+1);
	// RooDataSet set1("s1", "excited energy", fit_tree[1], RooArgSet(x1));
	// number of peaks
	size_t n1 = 10;
	// mean value
	const double mean1_value[16][3] = {
		{17.2, 17.0, 17.3},		// 17.2
		{0.0, 1.0, -1.0},		// 17.7
		{0.0, 1.0, -1.0},		// 18.5
		{0.0, 1.0, -1.0},		// 19.3
		{0.0, 1.0, -1.0},		// 20.1
		{0.0, 1.0, -1.0},		// 21.4
		{0.0, 1.0, -1.0},		// 22.4
		// {0.0, 1.0, -1.0},		// 23.0
		{23.6, 23.4, 23.9},		// 23.5
		{0.0, 1.0, -1.0},		// 24.5
		{0.0, 1.0, -1.0},		// 25.9
	};
	// g factor value
	const double gamma1_value[16][3] = {
		{0.05, 0.0, 1.0},		// 17.2
		{0.0, 1.0, -1.0},		// 17.7
		{0.0, 1.0, -1.0},		// 18.5
		{0.0, 1.0, -1.0},		// 19.3
		{0.0, 1.0, -1.0},		// 20.1
		{0.0, 1.0, -1.0},		// 21.4
		{0.0, 1.0, -1.0},		// 22.4
		// {0.0, 1.0, -1.0},		// 23.0
		{0.02, 0.0, 1.0},		// 23.5
		{0.0, 1.0, -1.0},		// 24.5
		{0.0, 1.0, -1.0},		// 25.9
	};
	// sigma value
	const double sigma1_value[16] = {
		0.209, 0.24, 0.28, 0.31,	// 17.2, 17.7, 18.5, 19.3
		0.34, 0.36, 0.38,	// 20.1, 21.4, 22.4, 23.0
		0.41, 0.42, 0.43,		// 23.5, 24.5, 25.9
	};
	// variable names
	const std::string name1[16] = {
		"172", "177", "185", "193",	// 17.3, 17.7, 18.5, 19.3
		"201", "214", "224", // 20.1, 21.4, 22.4, 23.0
		"235", "245", "259",		// 23.5, 24.5, 25.9
	};
	// threshold
	const double threshold1[16] = {
		// 17.2, 17.7, 18.5, 19.3
		threshold[1], threshold[1], threshold[1], threshold[1],
		// 20.1, 21.4, 22.4, 23.0
		threshold[2], threshold[2], threshold[2],
		// 23.5, 24.5, 25.9
		threshold[2], threshold[2], threshold[2]
	};
	// mean variable
	RooRealVar *mean1[16];
	// g factor variable
	RooRealVar *gamma1[16];
	// search from info0
	for (size_t i = 0; i < n1; ++i) {
		auto search = info_map[0].find(name1[i]);
		if (search == info_map[0].end()) {
			mean1[i] = nullptr;
			gamma1[i] = nullptr;
		} else {
			mean1[i] = search->second.mean;
			gamma1[i] = search->second.gamma;
		}
	}
	// sigma variable
	RooConstVar *sigma1[16];
	// P.D.F.
	RooAbsPdf *pdf1[16];
	// peak at other energy
	for (size_t i = 0; i < n1; ++i) {
		if (!mean1[i]) {
			mean1[i] = new RooRealVar(
				("mean"+name1[i]).c_str(), "mean",
				mean1_value[i][0], mean1_value[i][1], mean1_value[i][2]
			);
		}
		if (!gamma1[i]) {
			gamma1[i] = new RooRealVar(
				("gamma"+name1[i]).c_str(), "gamma",
				gamma1_value[i][0], gamma1_value[i][1], gamma1_value[i][2]
			);
		}
		sigma1[i] = new RooConstVar(
			("sigma"+name1[i]).c_str(), "sigma",
			sigma1_value[i]
		);
		pdf1[i] = new AsymmetricVoigtian(
			("voigt"+name1[i]).c_str(), "voigt",
			x, *mean1[i], *gamma1[i], *sigma1[i],
			threshold1[i], g_efficiency[1]
		);
		VoigtInfo info;
		info.mean = mean1[i];
		info.gamma = gamma1[i];
		info.sigma = sigma1[i];
		info.threshold = threshold1[i];
		info_map[1].insert(std::make_pair(name1[i], info));
	}
	// add fraction
	RooRealVar* frac1[16];
	// add result
	RooAbsPdf* model1[16];
	// add gaussian
	model1[0] = pdf1[0];
	for (size_t i = 1; i < n1; ++i) {
		frac1[i-1] = new RooRealVar(
			TString::Format("frac1_%ld", i), "frac",
			1.0-1.0/i, 0.0, 1.0
		);
		model1[i] = new RooAddPdf(
			TString::Format("model1_%ld", i), "model",
			RooArgList(*model1[i-1], *pdf1[i]),
			*frac1[i-1]
		);
	}

	// 10Be 6 MeV
	RooDataHist hist2("dh2", "excited energy", RooArgList(x), hist_ex+2);
	// RooDataSet set2("s2", "excited energy", fit_tree[2], RooArgSet(x));
	// number of peaks
	size_t n2 = 6;
	// mean value
	const double mean2_value[16][3] = {
		{0.0, 1.0, -1.0},		// 20.1
		{0.0, 1.0, -1.0},		// 21.4
		{0.0, 1.0, -1.0},		// 22.4
		// {0.0, 1.0, -1.0},		// 23.0
		{0.0, 1.0, -1.0},		// 23.5
								// 24.0
		{0.0, 1.0, -1.0},		// 24.5
		{0.0, 1.0, -1.0},		// 25.9
	};
	// g factor value
	const double gamma2_value[16][3] = {
		{0.0, 1.0, -1.0},		// 20.1
		{0.0, 1.0, -1.0},		// 21.4
		{0.0, 1.0, -1.0},		// 22.4
		// {0.0, 1.0, -1.0},		// 23.0
		{0.0, 1.0, -1.0},		// 23.5
		{0.0, 1.0, -1.0},		// 24.5
		{0.0, 1.0, -1.0},		// 25.9
	};
	// sigma value
	const double sigma2_value[16] = {
		0.3, 0.32, 0.34,	// 20.1, 21.4, 22.4, 23.0
		0.36, 0.37, 0.38,		// 23.5, 24.5, 25.9
	};
	// variable name
	const std::string name2[16] = {
		"201", "214", "224", // 20.1, 21.4, 22.4, 23.0
		"235", "245", "259",		// 23.5, 24.5, 25.9
	};
	// threshold
	const double threshold2[16] = {
		// 20.1, 21.4, 22.4, 23.0
		threshold[2], threshold[2], threshold[2],
		// 23.8, 24.5, 25.9
		threshold[2], threshold[2], threshold[2]
	};
	// mean variable
	RooRealVar *mean2[16];
	// g factor variable
	RooRealVar *gamma2[16];
	// search from info0 and info1
	for (size_t i = 0; i < n2; ++i) {
		auto search0 = info_map[0].find(name2[i]);
		if (search0 == info_map[0].end()) {
			auto search1 = info_map[1].find(name2[i]);
			if (search1 == info_map[1].end()) {
				mean2[i] = nullptr;
				gamma2[i] = nullptr;
			} else {
				mean2[i] = search1->second.mean;
				gamma2[i] = search1->second.gamma;
			}
		} else {
			mean2[i] = search0->second.mean;
			gamma2[i] = search0->second.gamma;
		}
	}
	// sigma variable
	RooConstVar *sigma2[16];
	// P.D.F.
	RooAbsPdf *pdf2[16];
	// peak at other energy
	for (size_t i = 0; i < n2; ++i) {
		if (!mean2[i]) {
			mean2[i] = new RooRealVar(
				("mean"+name2[i]).c_str(), "mean",
				mean2_value[i][0], mean2_value[i][1], mean2_value[i][2]
			);
		}
		if (!gamma2[i]) {
			gamma2[i] = new RooRealVar(
				("gamma"+name2[i]).c_str(), "gamma",
				gamma2_value[i][0], gamma2_value[i][1], gamma2_value[i][2]
			);
		}
		sigma2[i] = new RooConstVar(
			("sigma"+name2[i]).c_str(), "sigma",
			sigma2_value[i]
		);
		pdf2[i] = new AsymmetricVoigtian(
			("voigt"+name2[i]).c_str(), "voigt",
			x, *mean2[i], *gamma2[i], *sigma2[i],
			threshold2[i], g_efficiency[2]
		);
		VoigtInfo info;
		info.mean = mean2[i];
		info.gamma = gamma2[i];
		info.sigma = sigma2[i];
		info.threshold = threshold2[i];
		info_map[2].insert(std::make_pair(name2[i], info));
	}
	// add fraction
	RooRealVar* frac2[16];
	// add result
	RooAbsPdf* model2[16];
	// add pdf
	model2[0] = pdf2[0];
	for (size_t i = 1; i < n2; ++i) {
		frac2[i-1] = new RooRealVar(
			TString::Format("frac2_%ld", i), "frac",
			1.0-1.0/i, 0.0, 1.0
		);
		model2[i] = new RooAddPdf(
			TString::Format("model2_%ld", i), "model",
			RooArgList(*model2[i-1], *pdf2[i]),
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

	// print results
	std::cout << std::setw(6) << "state"
		<< std::setw(16) << "mean"
		<< std::setw(16) << "Gamma"
		<< std::setw(16) << "sigma" << "\n";
	const size_t print_n = 27;
	int print_state[print_n] = {
		0, 0, 0, 0,
		1, 0, 1, 0,
		1, 0, 1, 2,
		0, 1, 2, 0,
		1, 2, 0, 1,
		2, 0, 1, 2,
		0, 1, 2,
	};
	std::string print_name[print_n] = {
		"149", "154", "165", "177",
		"177", "185", "185", "193",
		"193", "201", "201", "201",
		"214", "214", "214", "224",
		"224", "224", "238", "235",
		"235", "245", "245", "245",
		"259", "259", "259"
	};
	for (size_t i = 0; i < print_n; ++i) {
		int state = print_state[i];
		auto search = info_map[state].find(print_name[i]);
		if (search == info_map[state].end()) {
			std::cerr << "Error: Could not find "
				<< print_name[i] << " in state " << state << "\n";
		}
		auto &info = search->second;
		double mean = info.mean->getVal();
		double g = info.gamma->getVal();
		double thres = info.threshold;
		double gamma = g * sqrt(mean - thres);
		std::cout << std::setw(6) << print_state[i]
			<< std::setw(16) << mean
			<< std::setw(16) << gamma
			<< std::setw(16) << info.sigma->getVal() << "\n";
	}

	// be quiet
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
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
			combined_hist, *pdf0[i]
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
			combined_hist, *pdf1[i]
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
			combined_hist, *pdf2[i]
		);
	}

	// save
	for (int i = 0; i < 3; ++i) hist_ex[i].Write();
	frame0->Write("rh0");
	frame1->Write("rh1");
	frame2->Write("rh2");
	// graph of resolution
	for (size_t i = 0; i < n0; ++i) {
		graph_resolution[0].AddPoint(
			mean0[i]->getVal(), sigma0_value[i]
		);
	}
	for (size_t i = 0; i < n1; ++i) {
		graph_resolution[1].AddPoint(
			mean1[i]->getVal(), sigma1_value[i]
		);
	}
	for (size_t i = 0; i < n2; ++i) {
		graph_resolution[2].AddPoint(
			mean2[i]->getVal(), sigma2_value[i]
		);
	}
	for (int i = 0; i < 3; ++i) {
		graph_resolution[i].Write(TString::Format("res%d", i));
	}
	// three lines graph
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

	TCanvas *c2 = new TCanvas("c2", "", 1920, 1080);
	test_frame->Draw();
	c2->SaveAs(
		TString::Format("%simage/test-asym-voigt.png", kGenerateDataPath)
	);
	c2->Write("tc");

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