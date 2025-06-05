#include <iostream>
#include <iomanip>
#include <map>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TLegend.h>
#include <RooRealVar.h>
#include <RooRealConstant.h>
#include <RooPlot.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooExtendPdf.h>
#include <RooGaussian.h>
#include <RooAddPdf.h>
#include <RooFitResult.h>
#include <RooCategory.h>
#include <RooSimultaneous.h>

#include "include/defs.h"
#include "include/spectrum/asymmetric_voigt.h"
#include "include/spectrum/background_poly.h"


using namespace ribll;

constexpr double threshold[3] = {12.01, 15.38, 18.19};
constexpr double resolution_fraction = 1.4;
constexpr double sigma_range = 0.0;
constexpr int calc_index = 0;

struct VoigtInfo {
	RooRealVar *mean;
	RooRealVar *gamma;
	RooConstVar *sigma;
	double threshold;
	RooRealVar *num;
};


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


struct SpectrumV3Event {
	int valid, ppac_flag, taf_flag;
	int be_state;
	double excited_energy;
};


int FillFromSpectrumV3(
	TH1F *hist,
	double *tree_ex,
	TTree **tree
) {
	// spectrum V3 file
	TString input_file_name = TString::Format(
		"%s%sC14-10Be-4He-2H-v3.root", kGenerateDataPath, kSpectrumDir
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
	SpectrumV3Event event;
	// setup input branches
	ipt->SetBranchAddress("valid", &event.valid);
	ipt->SetBranchAddress("be_state", &event.be_state);
	ipt->SetBranchAddress("excited_energy", &event.excited_energy);

	// loop
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		if (event.valid != 0) continue;

		if (event.be_state >= 0) {
			hist[event.be_state].Fill(event.excited_energy);
			tree_ex[event.be_state] = event.excited_energy;
			tree[event.be_state]->Fill();
		}
	}

	ipf.Close();
	return 0;
}


int main() {
	// efficiency file name
	TString efficiency_file_name = TString::Format(
		"%s%sefficiency-v2.root", kGenerateDataPath, kSimulateDir
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
		"%s%sroofit-v2.root", kGenerateDataPath, kSpectrumDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// histogram of all valid PPAC events
	TH1F hist_ex[3] = {
		TH1F("h0", "excited energy", 50, 12.06, 27.06),
		TH1F("h1", "excited energy", 50, 12.06, 27.06),
		TH1F("h2", "excited energy", 50, 12.06, 27.06)
	};

	// 分辨率
	TGraph graph_resolution[3];
	// 分支比
	TH1F relative_bar[3];
	relative_bar[0] = TH1F("hb0", "relative events", 11, 0, 11);
	relative_bar[0].SetBarOffset(0.55);
	relative_bar[0].SetBarWidth(0.3);
	relative_bar[0].SetFillColor(kBlue);
	relative_bar[1] = TH1F("hb1", "relative events", 11, 0, 11);
	relative_bar[1].SetBarOffset(0.85);
	relative_bar[1].SetBarWidth(0.3);
	relative_bar[1].SetFillColor(kGreen);
	relative_bar[2] = TH1F("hb2", "relative events", 11, 0, 11);
	relative_bar[2].SetBarOffset(1.15);
	relative_bar[2].SetBarWidth(0.3);
	relative_bar[2].SetFillColor(kRed);

	TTree *fit_tree[3];
	double fit_ex[3];
	for (int i = 0; i < 3; ++i) {
		fit_tree[i] = new TTree("tree", "tree for RooFit");
		fit_tree[i]->Branch("ex", fit_ex+i, "ex/D");
	}

	opf.cd();
	for (int i = 0; i < 3; ++i) {
		hist_ex[i].Write(TString::Format("ho%d", i));
	}
	if (FillFromSpectrumV3(hist_ex, fit_ex, fit_tree) != 0) {
		return -1;
	}
	opf.cd();
	for (int i = 0; i < 3; ++i) {
		hist_ex[i].Write(TString::Format("hoe%d", i));
	}

	// fitting function information
	std::map<std::string, VoigtInfo> info_map[3];

	// RooFit
	// get data
	RooRealVar x("ex", "", 12.06, 27.06);

	RooDataHist hist0("dh0", "", x, hist_ex);
	// RooDataSet set0("s0", "excited energy", fit_tree[0], RooArgSet(x));
	// construct P.D.F. information
	// 14.4
	info_map[0].insert(std::make_pair(
        "144",
        VoigtInfo {
            new RooRealVar("mean144", "mean", 14.4, 14.1, 14.5),
            new RooRealVar("gamma144", "gamma", 0.02, 1e-6, 0.2),
            new RooConstVar("sigma0_144", "sigma", 0.12*resolution_fraction),
            threshold[0],
			new RooRealVar("num0_144", "num", 4, 0.1, 10)
        }
    ));
	// 14.9
	info_map[0].insert(std::make_pair(
        "149",
        VoigtInfo {
            new RooRealVar("mean149", "mean", 14.9, 14.8, 15.1),
            new RooRealVar("gamma149", "gamma", 0.059, 0.02, 0.5),
            new RooConstVar("sigma0_149", "sigma", 0.14*resolution_fraction),
            threshold[0],
			new RooRealVar("num0_149", "num", 32, 10, 35)
        }
    ));
	// 15.6
	info_map[0].insert(std::make_pair(
        "156",
        VoigtInfo {
            new RooRealVar("mean154", "mean", 15.6, 15.5, 15.7),
            new RooRealVar("gamma154", "gamma", 0.095, 0.07, 0.2),
            new RooConstVar("sigma0_154", "sigma", 0.14*resolution_fraction),
            threshold[0],
			new RooRealVar("num0_154", "num", 30, 10, 45)
        }
    ));
    // 16.5
    info_map[0].insert(std::make_pair(
        "165",
        VoigtInfo {
            new RooRealVar("mean165", "mean", 16.3, 16.0, 16.6),
			new RooRealVar("gamma165", "gamma", 0.06, 0.05, 0.2),
            new RooConstVar("sigma0_165", "sigma", 0.16*resolution_fraction),
            threshold[0],
			new RooRealVar("num0_165", "num", 20, 2, 90)
		}
	));
	// 17.2
	info_map[0].insert(std::make_pair(
        "172",
        VoigtInfo {
            new RooRealVar("mean172", "mean", 17.2, 17.1, 17.4),
            new RooRealVar("gamma172", "gamma", 0.052, 0.04, 0.4),
            new RooConstVar("sigma0_172", "sigma", 0.18*resolution_fraction),
            threshold[1],
			new RooRealVar("num0_172", "num", 10, 0.01, 100)
        }
    ));
    // // 17.9
    // info_map[0].insert(std::make_pair(
    //     "179",
    //     VoigtInfo {
    //         new RooRealVar("mean179", "mean", 17.9, 17.6, 18.2),
    //         new RooRealVar("gamma179", "gamma", 0.04, 1e-6, 0.4),
	// 		new RooConstVar("sigma0_179", "sigma", 0.18*resolution_fraction),
    //         threshold[1],
	// 		new RooRealVar("num0_179", "num", 50, 20, 150)
	// 	}
	// ));
	// 18.2
	info_map[0].insert(std::make_pair(
        "182",
        VoigtInfo {
            new RooRealVar("mean182", "mean", 18.2, 18.0, 18.5),
            new RooRealVar("gamma182", "gamma", 0.02, 0.001, 0.4),
            new RooConstVar("sigma0_182", "sigma", 0.2*resolution_fraction),
            threshold[1],
			new RooRealVar("num0_182", "num", 50, 20, 150)
        }
    ));
	// 19.0
	info_map[0].insert(std::make_pair(
        "192",
        VoigtInfo {
            new RooRealVar("mean192", "mean", 19.2, 19.0, 19.4),
            new RooRealVar("gamma192", "gamma", 0.2, 0.05, 1.0),
            new RooConstVar("sigma0_192", "sigma", 0.21*resolution_fraction),
            threshold[1],
			new RooRealVar("num0_192", "num", 50, 1.0, 200)
        }
    ));
	// 19.8
	info_map[0].insert(std::make_pair(
        "198",
        VoigtInfo {
            new RooRealVar("mean198", "mean", 19.8, 19.6, 20.1),
            new RooRealVar("gamma198", "gamma", 1e-4, 1e-6, 0.5),
            new RooConstVar("sigma0_198", "sigma", 0.21*resolution_fraction),
            threshold[1],
			new RooRealVar("num0_198", "num", 20, 2, 110)
        }
    ));
    // 20.7
    info_map[0].insert(std::make_pair(
        "207",
        VoigtInfo {
            new RooRealVar("mean207", "mean", 20.6, 20.4, 20.9),
            new RooRealVar("gamma207", "gamma", 0.23, 0.02, 0.5),
			new RooConstVar("sigma0_207", "sigma", 0.23*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_207", "num", 1, 0.01, 100)
		}
	));
	// 21.6
	info_map[0].insert(std::make_pair(
		"216",
		VoigtInfo {
            new RooRealVar("mean216", "mean", 21.4, 21.3, 21.8),
            new RooRealVar("gamma216", "gamma", 0.12, 0.05, 0.5),
            new RooConstVar("sigma0_216", "sigma", 0.22*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_216", "num", 50, 0.01, 100)
		}
	));
	// 22.3
	info_map[0].insert(std::make_pair(
		"223",
		VoigtInfo {
            new RooRealVar("mean223", "mean", 22.3, 21.9, 22.5),
            new RooRealVar("gamma223", "gamma", 0.3, 0.01, 0.3),
            new RooConstVar("sigma0_223", "sigma", 0.24*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_223", "num", 10, 0.01, 40)
		}
	));
	// 22.8
	info_map[0].insert(std::make_pair(
		"228",
		VoigtInfo {
            new RooRealVar("mean228", "mean", 22.8, 22.7, 23.0),
            new RooRealVar("gamma228", "gamma", 0.05, 0.1, 0.3),
            new RooConstVar("sigma0_228", "sigma", 0.24*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_228", "num", 20, 0.01, 50)
		}
	));
	// 23.5
	info_map[0].insert(std::make_pair(
		"235",
		VoigtInfo {
            new RooRealVar("mean235", "mean", 23.5, 23.4, 23.6),
            new RooRealVar("gamma235", "gamma", 0.1, 0.01, 0.2),
            new RooConstVar("sigma0_235", "sigma", 0.25*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_235", "num", 20, 0.01, 50)
		}
	));
	// 24.2
	info_map[0].insert(std::make_pair(
		"245",
		VoigtInfo {
            new RooRealVar("mean245", "mean", 24.2, 23.9, 24.7),
            new RooRealVar("gamma245", "gamma", 0.01, 0.01, 0.2),
            new RooConstVar("sigma0_245", "sigma", 0.25*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_245", "num", 5, 0.01, 50)
		}
	));
	// // 24.8
	// info_map[0].insert(std::make_pair(
	// 	"248",
	// 	VoigtInfo {
    //         new RooRealVar("mean248", "mean", 24.9, 24.8, 25.5),
    //         new RooRealVar("gamma248", "gamma", 0.01, 0.01, 0.2),
    //         new RooConstVar("sigma0_248", "sigma", 0.25*resolution_fraction),
	// 		threshold[2],
	// 		new RooRealVar("num0_248", "num", 5, 0.01, 50)
	// 	}
	// ));
	// // 26.0
	// info_map[0].insert(std::make_pair(
	// 	"260",
	// 	VoigtInfo {
    //         new RooRealVar("mean260", "mean", 26.1, 25.6, 26.5),
    //         new RooRealVar("gamma260", "gamma", 0.5, 1e-6, 2.0),
    //         new RooConstVar("sigma0_260", "sigma", 0.26*resolution_fraction),
	// 		threshold[2],
	// 		new RooRealVar("num0_260", "num", 4, 0.01, 20)
	// 	}
	// ));


	// 10Be 3.5 MeV
	// RooRealVar x("ex", "ex1", 12.0, 27.0);
	RooDataHist hist1("dh1", "", x, hist_ex+1);
	// RooDataSet set1("s1", "excited energy", fit_tree[1], RooArgSet(x1));
	// construct P.D.F. information
	// 17.2
	info_map[1].insert(std::make_pair(
		"172",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_172", "sigma", 0.11*resolution_fraction),
            threshold[1],
			new RooRealVar("num1_172", "num", 30, 20, 150)
		}
	));
	// // 17.9
	// info_map[1].insert(std::make_pair(
	// 	"179",
	// 	VoigtInfo{
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma1_179", "sigma", 0.12*resolution_fraction),
    //         threshold[1],
	// 		new RooRealVar("num1_179", "num", 5, 0.01, 80)
	// 	}
	// ));
	// 18.2
	info_map[1].insert(std::make_pair(
		"182",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_182", "sigma", 0.13*resolution_fraction),
            threshold[1],
			new RooRealVar("num1_182", "num", 40, 20, 80)
		}
	));
	// 19.0
	info_map[1].insert(std::make_pair(
        "192",
        VoigtInfo {
			nullptr,
			nullptr,
            new RooConstVar("sigma1_192", "sigma", 0.15*resolution_fraction),
            threshold[1],
			new RooRealVar("num1_192", "num", 150, 50, 400)
        }
    ));
	// 19.8
	info_map[1].insert(std::make_pair(
        "198",
        VoigtInfo {
			nullptr,
			nullptr,
            new RooConstVar("sigma1_198", "sigma", 0.16*resolution_fraction),
            threshold[1],
			new RooRealVar("num1_198", "num", 100, 50, 150)
        }
    ));
	// 20.7
	info_map[1].insert(std::make_pair(
		"207",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_207", "sigma", 0.18*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_207", "num", 100, 50, 200)
		}
	));
	// 21.6
	info_map[1].insert(std::make_pair(
		"216",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_216", "sigma", 0.19*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_216", "num", 50, 20, 120)
		}
	));
	// 22.3
	info_map[1].insert(std::make_pair(
		"223",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_223", "sigma", 0.20*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_223", "num", 55, 50, 130)
		}
	));
	// 22.8
	info_map[1].insert(std::make_pair(
		"228",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_228", "sigma", 0.21*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_228", "num", 100, 10, 140)
		}
	));
	// 23.5
	info_map[1].insert(std::make_pair(
		"235",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_235", "sigma", 0.21*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_235", "num", 30, 0.1, 70)
		}
	));
	// 24.5
	info_map[1].insert(std::make_pair(
		"245",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_245", "sigma", 0.23*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_245", "num", 25, 10, 100)
		}
	));
	// // 24.8
	// info_map[1].insert(std::make_pair(
	// 	"248",
	// 	VoigtInfo{
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma1_248", "sigma", 0.23*resolution_fraction),
    //         threshold[2],
	// 		new RooRealVar("num1_248", "num", 25, 10, 100)
	// 	}
	// ));
	// // 26.0
	// info_map[1].insert(std::make_pair(
	// 	"260",
	// 	VoigtInfo{
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma1_260", "sigma", 0.25*resolution_fraction),
    //         threshold[2],
	// 		new RooRealVar("num1_260", "num", 100, 0.1, 200)
	// 	}
	// ));


	// 10Be 6 MeV
	RooDataHist hist2("dh2", "", RooArgList(x), hist_ex+2);
	// RooDataSet set2("s2", "excited energy", fit_tree[2], RooArgSet(x));
	// construct P.D.F. information
	// 19.8
	info_map[2].insert(std::make_pair(
		"198",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_198", "sigma", 0.11*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_198", "num", 30, 25, 60)
		}
	));
	// 20.7
	info_map[2].insert(std::make_pair(
		"207",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_207", "sigma", 0.11*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_207", "num", 80, 50, 150)
		}
	));
	// 21.6
	info_map[2].insert(std::make_pair(
		"216",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_216", "sigma", 0.14*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_216", "num", 30, 0.1, 150)
		}
	));
	// 22.3
	info_map[2].insert(std::make_pair(
		"223",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_223", "sigma", 0.15*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_223", "num", 200, 0.1, 400)
		}
	));
	// 22.8
	info_map[2].insert(std::make_pair(
		"228",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_228", "sigma", 0.17*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_228", "num", 10, 0.1, 150)
		}
	));
	// 23.5
	info_map[2].insert(std::make_pair(
		"235",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_235", "sigma", 0.18*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_235", "num", 10, 0.01, 100)
		}
	));
	// 24.5
	info_map[2].insert(std::make_pair(
		"245",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_245", "sigma", 0.18*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_245", "num", 60, 20, 100)
		}
	));
	// // 24.8
	// info_map[2].insert(std::make_pair(
	// 	"248",
	// 	VoigtInfo {
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma2_248", "sigma", 0.18*resolution_fraction),
	// 		threshold[2],
	// 		new RooRealVar("num2_248", "num", 60, 20, 100)
	// 	}
	// ));
	// // 26.0
	// info_map[2].insert(std::make_pair(
	// 	"260",
	// 	VoigtInfo {
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma2_260", "sigma", 0.23*resolution_fraction),
	// 		threshold[2],
	// 		new RooRealVar("num2_260", "num", 100, 0.1, 200)
	// 	}
	// ));


	// P.D.F
	std::vector<RooAbsPdf*> pdf0;
	// extended P.D.F.
	std::vector<RooExtendPdf*> extend_pdf0;
	for (auto &[key, value] : info_map[0]) {
		// check
		if (!value.mean || !value.gamma || !value.sigma) {
			std::cerr << "Error: Missing variable for state 0, energy "
				<< key << "\n";
			return -1;
		}
		// construct P.D.F.
		pdf0.push_back(new AsymmetricVoigtian(
			("voigt0_"+key).c_str(), "voigt",
			x,
			*value.mean, *value.gamma, *value.sigma,
			value.threshold,
			g_efficiency[0]
		));
		// construct extended P.D.F.
		extend_pdf0.push_back(new RooExtendPdf(
			("ext_voigt0_"+key).c_str(), "extend",
            *(pdf0.back()), *value.num
        ));
	}

	// arguments list for model0
	RooArgList model0_list;
	// construct extended P.D.F.
	for (size_t i = 0; i < pdf0.size(); ++i) {
		model0_list.add(*extend_pdf0[i]);
	}
	// background
	RooRealVar bkg0_root0("bkg0_root0", "root0", 14.0, 14.0, 16.0);
	RooRealVar bkg0_root1("bkg0_root1", "root1", 26.0, 25.0, 27.0);
	BackgroundPoly bkg0("bkg0", "bkg", x, bkg0_root0, bkg0_root1);
	// background events' number
	RooRealVar num_bkg0("num_bkg0", "num", 10, 0, 250);
	RooExtendPdf ext_bkg0("ext_bkg0", "extend", bkg0, num_bkg0);
	model0_list.add(ext_bkg0);
	// model0
	RooAddPdf model0("model0", "model", model0_list);


	// P.D.F.
	std::vector<RooAbsPdf*> pdf1;
	// extend P.D.F.
	std::vector<RooExtendPdf*> extend_pdf1;
	for (auto &[key, value] : info_map[1]) {
		// get same mean and gamma values
		auto search = info_map[0].find(key);
		if (search != info_map[0].end()) {
			value.mean = search->second.mean;
			value.gamma = search->second.gamma;
		}
		// check
		if (!value.mean || !value.gamma || !value.sigma) {
			std::cerr << "Error: Missing variable for state 1, energy "
				<< key << "\n";
			return -1;
		}
		// construct P.D.F.
		pdf1.push_back(new AsymmetricVoigtian(
			("voigt1_"+key).c_str(), "voigt",
			x,
			*value.mean, *value.gamma, *value.sigma,
            value.threshold,
			g_efficiency[1]
		));
		// construct extend P.D.F.
		extend_pdf1.push_back(new RooExtendPdf(
			("ext_voigt1_"+key).c_str(), "extend",
			*(pdf1.back()), *value.num
		));
	}
	// argument list for model1
	RooArgList model1_list;
	// construct extend P.D.F.
	for (size_t i = 0; i < pdf1.size(); ++i) {
		model1_list.add(*extend_pdf1[i]);
	}
	RooRealVar tail1_mean("tail1_mean", "mean", 26.0, 25.5, 26.5);
	RooRealVar tail1_sigma("tail1_sigma", "sigma", 0.5, 0.001, 1.0);
	RooGaussian tail1_gaus("tail1_gaus", "gaus", x, tail1_mean, tail1_sigma);
	RooRealVar tail1_num("num_tail1", "num", 20, 0.1, 100);
	RooExtendPdf ext_tail1("ext_tail1", "extend", tail1_gaus, tail1_num);
	model1_list.add(ext_tail1);
	// background
	RooRealVar bkg1_root0("bkg1_root0", "root0", 17.0, 16.0, 18.0);
	RooRealVar bkg1_root1("bkg1_root1", "root1", 34.0, 33.0, 35.0);
	BackgroundPoly bkg1("bkg1", "bkg", x, bkg1_root0, bkg1_root1);
	// background events' number
	RooRealVar num_bkg1("num_bkg1", "num", 10, 0.1, 250);
	RooExtendPdf ext_bkg1("ext_bkg1", "extend", bkg1, num_bkg1);
	model1_list.add(ext_bkg1);
	// model1
	RooAddPdf model1("model1", "model", model1_list);


	// P.D.F.
	std::vector<RooAbsPdf*> pdf2;
	// extended P.D.F.
	std::vector<RooExtendPdf*> extend_pdf2;
	for (auto &[key, value] : info_map[2]) {
		// get same mean and gamma values
		auto search = info_map[0].find(key);
		if (search != info_map[0].end()) {
			value.mean = search->second.mean;
			value.gamma = search->second.gamma;
		} else {
			auto search1 = info_map[1].find(key);
			if (search1 != info_map[1].end()) {
				value.mean = search1->second.mean;
				value.gamma = search1->second.gamma;
			}
		}
		// check
		if (!value.mean || !value.gamma || !value.sigma) {
			std::cerr << "Error: Missing variable for state 2, energy "
				<< key << "\n";
			return -1;
		}
		// construct P.D.F.
		pdf2.push_back(new AsymmetricVoigtian(
			("voigt2_"+key).c_str(), "voigt",
			x,
			*value.mean, *value.gamma, *value.sigma,
            value.threshold,
			g_efficiency[2]
		));
		// construct extend P.D.F.
		extend_pdf2.push_back(new RooExtendPdf(
            ("ext_voigt2_"+key).c_str(), "extend",
            *(pdf2.back()), *value.num
        ));
	}
	// arguments list for model0
	RooArgList model2_list;
	// construct extended P.D.F.
	for (size_t i = 0; i < pdf2.size(); ++i) {
		model2_list.add(*extend_pdf2[i]);
	}
	RooRealVar tail2_mean("tail2_mean", "mean", 26.5, 25.8, 27.0);
	RooRealVar tail2_sigma("tail2_sigma", "sigma", 0.5, 0.01, 1.0);
	RooGaussian tail2_gaus("tail2_gaus", "gaus", x, tail2_mean, tail2_sigma);
	RooRealVar tail2_num("num_tail2", "num", 20, 0.1, 100);
	RooExtendPdf ext_tail2("ext_tail2", "extend", tail2_gaus, tail2_num);
	model2_list.add(ext_tail2);
	// background
	RooRealVar bkg2_root0("bkg2_root0", "root0", 20.0, 19.5, 21.0);
	RooRealVar bkg2_root1("bkg2_root1", "root1", 34.0, 33.0, 35.0);
	BackgroundPoly bkg2("bkg2", "bkg", x, bkg2_root0, bkg2_root1);
	// background events' number
	RooRealVar num_bkg2("num_bkg2", "num", 10, 0.1, 30);
	RooExtendPdf ext_bkg2("ext_bkg2", "extend", bkg2, num_bkg2);
	model2_list.add(ext_bkg2);
	// model0
	RooAddPdf model2("model2", "model", model2_list);


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
		{
			{"0", &model0},
			{"1", &model1},
			{"2", &model2}
		},
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
	std::cout << "NLL: " << fit_result->minNll() << "\n";

	// show bar
	std::vector<std::string> bar_names = {
		"172", "182", "192",
		"198", "207", "216", "223",
		"228", "235", "245",
	};
	for (size_t i = 0; i < bar_names.size(); ++i) {
		for (int state = 0; state < 3; ++state) {
			auto search = info_map[state].find(bar_names[i]);
			if (search == info_map[state].end()) {
				std::cerr << "Warning: Could not find "
					<< bar_names[i] << " in state " << state << "\n";
				continue;
			}
			double mean = search->second.mean->getVal();
			double event = search->second.num->getVal();
			double relative_event = event / g_efficiency[state]->Eval(mean);
			relative_bar[state].SetBinContent(i+1, relative_event);
		}
	}

	// print results
	std::cout << std::setw(6) << "state"
		<< std::setw(16) << "mean"
		<< std::setw(16) << "Gamma"
		<< std::setw(16) << "sigma"
		<< std::setw(16) << "events"
		<< std::setw(16) << "relative events" << "\n";
	std::vector<std::pair<int, std::string>> print_info = {
		{0, "144"},
		{0, "149"},
		{0, "156"},
		{0, "165"},
		{0, "172"},
		{1, "172"},
		// {0, "179"},
		// {1, "179"},
		{0, "182"},
		{1, "182"},
		{0, "192"},
		{1, "192"},
		{0, "198"},
		{1, "198"},
		{2, "198"},
		{0, "207"},
		{1, "207"},
		{2, "207"},
		{0, "216"},
		{1, "216"},
		{2, "216"},
		{0, "223"},
		{1, "223"},
		{2, "223"},
		{0, "228"},
		{1, "228"},
		{2, "228"},
		{0, "235"},
        {1, "235"},
        {2, "235"},
		{0, "245"},
		{1, "245"},
		{2, "245"},
		// {0, "248"},
		// {1, "248"},
		// {2, "248"},
	};

	for (size_t i = 0; i < print_info.size(); ++i) {
		int state = print_info[i].first;
		auto search = info_map[state].find(print_info[i].second);
		if (search == info_map[state].end()) {
			std::cerr << "Error: Could not find "
				<< print_info[i].second << " in state " << state << "\n";
			continue;
		}
		auto &info = search->second;
		double mean = info.mean->getVal();
		double g = info.gamma->getVal();
		double thres = info.threshold;
		double gamma = g * sqrt(mean - thres);
		double event = info.num->getVal();
		double relative_event = event / g_efficiency[state]->Eval(mean);
		std::cout << std::setw(6) << print_info[i].first
			<< std::setw(16) << mean
			<< std::setw(16) << gamma
			<< std::setw(16) << info.sigma->getVal()
			<< std::setw(16) << event
			<< std::setw(16) << relative_event << "\n";
	}
	// be quiet
	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
	// plot
	// frame 0
	RooPlot *frame0 = x.frame();
	combined_hist.plotOn(
		frame0,
		RooFit::Cut("state==state::0"),
		RooFit::DrawOption("0"),
		// RooFit::FillColor(kCyan),
		RooFit::DataError(RooAbsData::None)
	);
	simultaneous_model.plotOn(
		frame0,
		RooFit::Slice(state, "0"),
		RooFit::ProjWData(state, combined_hist),
		RooFit::LineColor(kRed)
	);
	for (size_t i = 0; i < pdf0.size(); ++i) {
		PlotPdf(
			simultaneous_model, frame0,
			state, "0",
			combined_hist, *pdf0[i]
		);
	}
	simultaneous_model.plotOn(
		frame0,
		RooFit::Slice(state, "0"),
		RooFit::ProjWData(state, combined_hist),
		RooFit::LineColor(kGreen),
		RooFit::Components(bkg0)
	);

	// frame 1
	RooPlot *frame1 = x.frame();
	combined_hist.plotOn(
		frame1,
		RooFit::Cut("state==state::1"),
		RooFit::DrawOption("0"),
		// RooFit::FillColor(kCyan),
		RooFit::DataError(RooAbsData::None)
	);
	simultaneous_model.plotOn(
		frame1,
		RooFit::Slice(state, "1"),
		RooFit::ProjWData(state, combined_hist),
		RooFit::LineColor(kRed)
	);
	for (size_t i = 0; i < pdf1.size(); ++i) {
		PlotPdf(
			simultaneous_model, frame1,
			state, "1",
			combined_hist, *pdf1[i]
		);
	}
	simultaneous_model.plotOn(
		frame1,
		RooFit::Slice(state, "1"),
		RooFit::ProjWData(state, combined_hist),
		RooFit::LineColor(kBlue),
		RooFit::Components(tail1_gaus)
	);
	simultaneous_model.plotOn(
		frame1,
		RooFit::Slice(state, "1"),
		RooFit::ProjWData(state, combined_hist),
		RooFit::LineColor(kGreen),
		RooFit::Components(bkg1)
	);

	// frame 2
	RooPlot *frame2 = x.frame();
	combined_hist.plotOn(
		frame2,
		RooFit::Cut("state==state::2"),
		RooFit::DrawOption("0"),
		// RooFit::FillColor(kCyan),
		RooFit::DataError(RooAbsData::None)
	);
	simultaneous_model.plotOn(
		frame2,
		RooFit::Slice(state, "2"),
		RooFit::ProjWData(state, combined_hist),
		RooFit::LineColor(kRed)
	);
	for (size_t i = 0; i < pdf2.size(); ++i) {
		PlotPdf(
			simultaneous_model, frame2,
			state, "2",
			combined_hist, *pdf2[i]
		);
	}
	simultaneous_model.plotOn(
		frame2,
		RooFit::Slice(state, "2"),
		RooFit::ProjWData(state, combined_hist),
		RooFit::LineColor(kBlue),
		RooFit::Components(tail2_gaus)
	);
	simultaneous_model.plotOn(
		frame2,
		RooFit::Slice(state, "2"),
		RooFit::ProjWData(state, combined_hist),
		RooFit::LineColor(kGreen),
		RooFit::Components(bkg2)
	);

	// printable efficiency
	TGraph print_efficiency[3];
	const double efficiency_factor[3] = {
		350.0, 600.0, 400.0
	};
	for (int i = 0; i < 3; ++i) {
		print_efficiency[i].SetLineWidth(2);
		print_efficiency[i].SetLineStyle(7);
		for (int j = 0; j < g_efficiency[i]->GetN(); ++j) {
			print_efficiency[i].AddPoint(
				g_efficiency[i]->GetPointX(j),
				g_efficiency[i]->GetPointY(j) * efficiency_factor[i]
			);
		}
	}

	// save
	opf.cd();
	for (int i = 0; i < 3; ++i) hist_ex[i].Write();
	frame0->Write("rh0");
	frame1->Write("rh1");
	frame2->Write("rh2");
	// // graph of resolution
	// for (size_t i = 0; i < n0-1; ++i) {
	// 	graph_resolution[0].AddPoint(
	// 		mean0[i]->getVal(), sigma0_value[i]
	// 	);
	// }
	for (size_t i = 0; i < 3; ++i) relative_bar[i].Write();
	for (const auto &[key, value] : info_map[0]) {
		graph_resolution[0].AddPoint(
            value.mean->getVal(), value.sigma->getVal()
        );
	}
	for (const auto &[key, value] : info_map[1]) {
		graph_resolution[1].AddPoint(
			value.mean->getVal(), value.sigma->getVal()
		);
	}
	for (const auto &[key, value] : info_map[2]) {
		graph_resolution[2].AddPoint(
			value.mean->getVal(), value.sigma->getVal()
		);
	}
	for (int i = 0; i < 3; ++i) {
		graph_resolution[i].Write(TString::Format("res%d", i));
	}
	// three lines graph
	gStyle->SetPadBorderMode(0);
	TCanvas *c1 = new TCanvas("c1", "", 1520, 885);
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
	// gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	frame0->GetYaxis()->SetLabelSize(0.15);
	frame0->SetTitle("");
	frame0->Draw();
	hist_ex[0].GetYaxis()->SetLabelSize(0.15);
	hist_ex[0].SetLineColor(kBlack);
	hist_ex[0].Draw("same");
	std::map<std::string, double> line0_y;
	line0_y.insert(std::make_pair("144", 5.0));
	line0_y.insert(std::make_pair("149", 12.0));
	line0_y.insert(std::make_pair("156", 20.0));
	line0_y.insert(std::make_pair("165", 12.0));
	line0_y.insert(std::make_pair("172", 20.0));
	line0_y.insert(std::make_pair("182", 30.0));
	line0_y.insert(std::make_pair("192", 15.0));
	line0_y.insert(std::make_pair("198", 24.0));
	line0_y.insert(std::make_pair("207", 10.0));
	line0_y.insert(std::make_pair("216", 11.0));
	line0_y.insert(std::make_pair("223", 5.0));
	line0_y.insert(std::make_pair("228", 10.0));
	line0_y.insert(std::make_pair("235", 6.0));
	line0_y.insert(std::make_pair("245", 5.0));
	for (const auto &[key, value] : info_map[0]) {
		// check y
		auto search = line0_y.find(key);
		double ymax = search == line0_y.end() ? 35.0  : search->second;
		TLine *line = new TLine(
			value.mean->getVal(), 0.0, value.mean->getVal(), ymax
		);
		line->SetLineStyle(2);
		line->Draw("same");
	}
	print_efficiency[0].Draw("same");

	pads[1]->cd();
	frame1->SetTitle("");
	frame1->Draw();
	frame1->GetYaxis()->SetLabelSize(0.15);
	hist_ex[1].SetLineColor(kBlack);
	hist_ex[1].GetYaxis()->SetLabelSize(0.15);
	hist_ex[1].Draw("same");
	print_efficiency[1].Draw("same");
	for (const auto &[key, value] : info_map[1]) {
		TLine *line = new TLine(
			value.mean->getVal(), 0.0, value.mean->getVal(), 100.0
		);
		line->SetLineStyle(2);
		line->Draw("same");
	}

	pads[2]->cd();
	frame2->SetTitle("");
	frame2->Draw();
	frame2->GetXaxis()->SetLabelSize(0.12);
	frame2->GetYaxis()->SetLabelSize(0.15);
	hist_ex[2].SetLineColor(kBlack);
	hist_ex[2].GetXaxis()->SetLabelSize(0.12);
	hist_ex[2].GetYaxis()->SetLabelSize(0.15);
	hist_ex[2].Draw("same");
	print_efficiency[2].Draw("same");
	for (const auto &[key, value] : info_map[2]) {
		TLine *line = new TLine(
			value.mean->getVal(), 0.0, value.mean->getVal(), 80.0
		);
		line->SetLineStyle(2);
		line->Draw("same");
	}

	c1->SaveAs(
		TString::Format("%s%sroofit-result-v2.png", kGenerateDataPath, kImageDir)
	);
	c1->Write("c1");

	TCanvas *c2 = new TCanvas("c2", "", 1920, 1080);
	gStyle->SetOptStat(0);
	TAxis *xaxis = relative_bar[1].GetXaxis();
	for (size_t i = 0; i < bar_names.size(); i++) {
		// get same mean and gamma values
		auto search = info_map[0].find(bar_names[i]);
		if (search == info_map[0].end()) continue;
		double mean = search->second.mean->getVal();
		xaxis->SetNdivisions(bar_names.size()+1);
		xaxis->ChangeLabel(
			i+2, -1, -1, -1, -1, -1,
			TString::Format("%.1lf", mean)
		);
	}
	xaxis->ChangeLabel(1, -1, 0, -1, -1, -1, "");
	xaxis->ChangeLabel(bar_names.size()+2, -1, 0, -1, -1, -1, "");
	relative_bar[1].Draw("bar");
	relative_bar[0].Draw("bar, same");
	relative_bar[2].Draw("bar, same");
	TLegend *legend = new TLegend(0.72, 0.75, 0.88, 0.85);
	legend->AddEntry(relative_bar+0, "ground state", "f");
	legend->AddEntry(relative_bar+1, "2+ state", "f");
	legend->AddEntry(relative_bar+2, "~6MeV state", "f");
	legend->Draw();
	c2->SaveAs(
		TString::Format("%s%srelative-events-v2.png", kGenerateDataPath, kImageDir)
	);
	c2->Write("barc");

	opf.Close();

	return 0;
}