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
constexpr double resolution_fraction = 1.2;
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
	int ppac_flag, taf_flag, target_flag, bind, hole, straight[2];
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
	ipt->SetBranchAddress("straight", straight);
	ipt->SetBranchAddress("q", q);
	ipt->SetBranchAddress("be_state", be_state);
	ipt->SetBranchAddress("excited_energy_target", excited_energy_target);
	ipt->SetBranchAddress("stateless_excited_energy", stateless_excited_energy);

	// efficiency file name
	TString efficiency_file_name = TString::Format(
		"%s%sefficiency-0002.root", kGenerateDataPath, kSimulateDir
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
	hist_ex[0] = TH1F("h0", "excited energy", 50, 12.10, 27.10);
	hist_ex[1] = TH1F("h1", "excited energy", 50, 12.10, 27.10);
	hist_ex[2] = TH1F("h2", "excited energy", 50, 12.10, 27.10);
	TGraph graph_resolution[3];
	TH1F relative_bar[3];
	relative_bar[0] = TH1F("hb0", "relative events", 13, 0, 13);
	relative_bar[0].SetBarOffset(0.55);
	relative_bar[0].SetBarWidth(0.3);
	relative_bar[0].SetFillColor(kBlue);
	relative_bar[1] = TH1F("hb1", "relative events", 13, 0, 13);
	relative_bar[1].SetBarOffset(0.85);
	relative_bar[1].SetBarWidth(0.3);
	relative_bar[1].SetFillColor(kGreen);
	relative_bar[2] = TH1F("hb2", "relative events", 13, 0, 13);
	relative_bar[2].SetBarOffset(1.15);
	relative_bar[2].SetBarWidth(0.3);
	relative_bar[2].SetFillColor(kRed);

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
		// straight PID
		// if (straight[0] != 3 || straight[1] != 3) continue;
		if (straight[0] != 3) continue;
		// if (
		// 	q[0] < -10 && q[0] > -14
		// 	&& q[2] < -10.15 && q[2] > -14.15
		// 	&& q[3] < -10 && q[3] > -14
		// ) {
		// 	hist_ex[0].Fill(stateless_excited_energy[0][calc_index]);
		// 	fit_ex[0] = stateless_excited_energy[0][calc_index];
		// 	fit_tree[0]->Fill();
		// }
		// if (
		// 	q[0] < -14.5 && q[0] > -16
		// 	&& q[2] < -14.65 && q[2] > -16.15
		// 	&& q[3] < -14.5 && q[3] > -16
		// ) {
		// 	hist_ex[1].Fill(stateless_excited_energy[1][calc_index]);
		// 	fit_ex[1] = stateless_excited_energy[1][calc_index];
		// 	fit_tree[1]->Fill();
		// }
		// if (
		// 	q[0] < -17.5 && q[0] > -20.0
		// 	&& q[2] < -17.65 && q[2] > -20.15
		// 	&& q[3] < -17.5 && q[3] > -20.0
		// ) {
		// 	hist_ex[2].Fill(stateless_excited_energy[2][calc_index]);
		// 	fit_ex[2] = stateless_excited_energy[2][calc_index];
		// 	fit_tree[2]->Fill();
		// }
		if (q[calc_index] < -11 && q[calc_index] > -13) {
			hist_ex[0].Fill(stateless_excited_energy[0][calc_index]);
			fit_ex[0] = stateless_excited_energy[0][calc_index];
			fit_tree[0]->Fill();
		} else if (q[calc_index] < -14.5 && q[calc_index] > -16) {
			hist_ex[1].Fill(stateless_excited_energy[1][calc_index]);
			fit_ex[1] = stateless_excited_energy[1][calc_index];
			fit_tree[1]->Fill();
		} else if (q[calc_index] < -17 && q[calc_index] > -20) {
			hist_ex[2].Fill(stateless_excited_energy[2][calc_index]);
			fit_ex[2] = stateless_excited_energy[2][calc_index];
			fit_tree[2]->Fill();
		}
	}

	for (int i = 0; i < 3; ++i) hist_ex[i].Write(TString::Format("ho%d", i));


	// test asymmetric voigt
	RooRealVar test_x("tx", "test x", 12.0, 27.0);
	RooRealVar test_mean("tmean", "test  mean", 16.0);
	RooRealVar test_g("tg", "test g", 0.2);
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
	RooRealVar x("ex", "excited energy", 12.10, 27.10);

	RooDataHist hist0("dh0", "excited energy", x, hist_ex);
	// RooDataSet set0("s0", "excited energy", fit_tree[0], RooArgSet(x));
	// construct P.D.F. information
	// 13.6
	info_map[0].insert(std::make_pair(
		"136",
		VoigtInfo {
			new RooRealVar("mean136", "mean", 13.65, 13.6, 13.7),
			new RooRealVar("gamma136", "gamma", 0.05, 0.0, 0.2),
			new RooConstVar("sigma0_136", "sigma", 0.1*resolution_fraction),
			threshold[0],
			new RooRealVar("num0_136", "num", 1, 0.01, 10)
		}
	));
	// 14.9
	info_map[0].insert(std::make_pair(
        "149",
        VoigtInfo {
            new RooRealVar("mean149", "mean", 14.9, 14.5, 15.2),
            new RooRealVar("gamma149", "gamma", 0.06, 0.0, 0.2),
            new RooConstVar("sigma0_149", "sigma", 0.14*resolution_fraction),
            threshold[0],
			new RooRealVar("num0_149", "num", 10, 0.01, 50)
        }
    ));
	// 15.4
	info_map[0].insert(std::make_pair(
        "154",
        VoigtInfo {
            new RooRealVar("mean154", "mean", 15.4, 15.3, 15.8),
            new RooRealVar("gamma154", "gamma", 0.1, 0.0, 0.2),
            new RooConstVar("sigma0_154", "sigma", 0.14*resolution_fraction),
            threshold[0],
			new RooRealVar("num0_154", "num", 15, 0.01, 50)
        }
    ));
    // 16.5
    info_map[0].insert(std::make_pair(
        "165",
        VoigtInfo {
            new RooRealVar("mean165", "mean", 16.5, 16.0, 17.0),
			new RooRealVar("gamma165", "gamma", 0.06, 0.0, 0.2),
            new RooConstVar("sigma0_165", "sigma", 0.16*resolution_fraction),
            threshold[0],
			new RooRealVar("num0_165", "num", 15, 0.01, 50)
		}
	));
	// 17.2
	info_map[0].insert(std::make_pair(
        "172",
        VoigtInfo {
            new RooRealVar("mean172", "mean", 17.2, 17.0, 17.5),
            new RooRealVar("gamma172", "gamma", 0.1, 0.0, 0.4),
            new RooConstVar("sigma0_172", "sigma", 0.18*resolution_fraction),
            threshold[1],
			new RooRealVar("num0_172", "num", 20, 0.01, 50)
        }
    ));
    // // 17.7
    // info_map[0].insert(std::make_pair(
    //     "177",
    //     VoigtInfo {
    //         new RooRealVar("mean177", "mean", 17.7, 17.4, 17.8),
    //         new RooRealVar("gamma177", "gamma", 0.05, 0.0, 0.2),
	// 		new RooConstVar("sigma0_177", "sigma", 0.19*resolution_fraction),
    //         threshold[1],
	// 		new RooRealVar("num0_136", "num", 1, 0.01, 10)
	// 	}
	// ));
	// 18.5
	info_map[0].insert(std::make_pair(
        "185",
        VoigtInfo {
            new RooRealVar("mean185", "mean", 18.5, 18.0, 18.6),
            new RooRealVar("gamma185", "gamma", 0.05, 0.0, 0.2),
            new RooConstVar("sigma0_185", "sigma", 0.2*resolution_fraction),
            threshold[1],
			new RooRealVar("num0_185", "num", 30, 0.01, 100)
        }
    ));
	// 19.0
	info_map[0].insert(std::make_pair(
        "190",
        VoigtInfo {
            new RooRealVar("mean190", "mean", 19.0, 18.5, 19.2),
            new RooRealVar("gamma190", "gamma", 0.05, 0.0, 0.4),
            new RooConstVar("sigma0_190", "sigma", 0.2*resolution_fraction),
            threshold[1],
			new RooRealVar("num0_190", "num", 1, 0.01, 10)
        }
    ));
    // 19.3
    info_map[0].insert(std::make_pair(
        "193",
        VoigtInfo {
            new RooRealVar("mean193", "mean", 19.3, 19.1, 19.6),
            new RooRealVar("gamma193", "gamma", 0.05, 0.0, 0.2),
			new RooConstVar("sigma0_193", "sigma", 0.21*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_193", "num", 30, 0.01, 100)
		}
	));
	// 20.1
	info_map[0].insert(std::make_pair(
        "201",
        VoigtInfo {
            new RooRealVar("mean201", "mean", 20.1, 19.7, 20.2),
            new RooRealVar("gamma201", "gamma", 0.02, 0.0, 0.4),
            new RooConstVar("sigma0_201", "sigma", 0.22*resolution_fraction),
            threshold[2],
			new RooRealVar("num0_201", "num", 20, 0.01, 50)
        }
    ));
    // 20.5
    info_map[0].insert(std::make_pair(
        "205",
        VoigtInfo {
            new RooRealVar("mean205", "mean", 20.5, 20.3, 20.7),
            new RooRealVar("gamma205", "gamma", 0.02, 0.0, 0.4),
			new RooConstVar("sigma0_205", "sigma", 0.23*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_205", "num", 30, 0.01, 50)
		}
	));
	// 20.9
	info_map[0].insert(std::make_pair(
		"209",
		VoigtInfo {
			new RooRealVar("mean209", "mean", 20.9, 20.7, 21.3),
            new RooRealVar("gamma209", "gamma", 0.05, 0.0, 0.2),
            new RooConstVar("sigma0_209", "sigma", 0.23*resolution_fraction),
            threshold[2],
			new RooRealVar("num0_209", "num", 15, 0.01, 50)
		}
	));
	// 21.6
	info_map[0].insert(std::make_pair(
		"216",
		VoigtInfo {
            new RooRealVar("mean216", "mean", 21.6, 21.3, 21.8),
            new RooRealVar("gamma216", "gamma", 0.05, 0.0, 0.2),
            new RooConstVar("sigma0_216", "sigma", 0.24*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_216", "num", 15, 0.01, 50)
		}
	));
	// 22.3
	info_map[0].insert(std::make_pair(
		"223",
		VoigtInfo {
            new RooRealVar("mean223", "mean", 22.3, 22.0, 22.6),
            new RooRealVar("gamma223", "gamma", 0.05, 0.0, 0.2),
            new RooConstVar("sigma0_223", "sigma", 0.25*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_223", "num", 5, 0.01, 20)
		}
	));
	// 22.8
	info_map[0].insert(std::make_pair(
		"228",
		VoigtInfo {
            new RooRealVar("mean228", "mean", 22.8, 22.7, 23.1),
            new RooRealVar("gamma228", "gamma", 0.05, 0.0, 0.2),
            new RooConstVar("sigma0_228", "sigma", 0.26*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_228", "num", 20, 0.01, 50)
		}
	));
	// 23.2
	info_map[0].insert(std::make_pair(
		"232",
		VoigtInfo {
            new RooRealVar("mean232", "mean", 23.2, 23.0, 23.5),
            new RooRealVar("gamma232", "gamma", 0.02, 0.0, 0.2),
            new RooConstVar("sigma0_232", "sigma", 0.26*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_232", "num", 5, 0.01, 50)
		}
	));
	// 23.6
	info_map[0].insert(std::make_pair(
		"236",
		VoigtInfo {
            new RooRealVar("mean236", "mean", 23.6, 23.3, 24.0),
            new RooRealVar("gamma236", "gamma", 0.02, 0.0, 0.2),
            new RooConstVar("sigma0_236", "sigma", 0.27*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_236", "num", 10, 0.01, 50)
		}
	));
	// 24.5
	info_map[0].insert(std::make_pair(
		"245",
		VoigtInfo {
            new RooRealVar("mean245", "mean", 24.5, 24.2, 24.7),
            new RooRealVar("gamma245", "gamma", 0.02, 0.0, 0.2),
            new RooConstVar("sigma0_245", "sigma", 0.28*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_245", "num", 1, 0.01, 10)
		}
	));
	// 25.4
	info_map[0].insert(std::make_pair(
		"254",
		VoigtInfo {
            new RooRealVar("mean254", "mean", 25.4, 25.2, 25.7),
            new RooRealVar("gamma254", "gamma", 0.02, 0.0, 1.0),
            new RooConstVar("sigma0_254", "sigma", 0.28*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_254", "num", 1, 0.01, 10)
		}
	));
	// 26.0
	info_map[0].insert(std::make_pair(
		"260",
		VoigtInfo {
            new RooRealVar("mean260", "mean", 26.2, 26.0, 26.4),
            new RooRealVar("gamma260", "gamma", 0.02, 0.0, 1.0),
            new RooConstVar("sigma0_260", "sigma", 0.31*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_260", "num", 1, 0.01, 10)
		}
	));

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

	// // number of events
	// std::vector<RooRealVar*> num_events0;
	// // extended P.D.F.
	// std::vector<RooExtendPdf*> extend_pdf0;
	// arguments list for model0
	RooArgList model0_list;
	// construct extended P.D.F.
	for (size_t i = 0; i < pdf0.size(); ++i) {
		// num_events0.push_back(new RooRealVar(
		// 	TString::Format("num0_%ld", i), "num",
		// 	10, 0, 1000
		// ));
		// extend_pdf0.push_back(new RooExtendPdf(
        //     TString::Format("extend0_%ld", i), "extend",
        //     *pdf0[i], *num_events0[i]
        // ));
		model0_list.add(*extend_pdf0[i]);
	}
	// background
	RooRealVar bkg0_p1("bkg0_p1", "p1", 40, 30, 60);
	RooRealVar bkg0_p0("bkg0_p0", "p0", -300, -1000, 0.0);
	BackgroundPoly bkg0("bkg0", "bkg", x, bkg0_p1, bkg0_p0);
	// background events' number
	RooRealVar num_bkg0("num_bkg0", "num", 10, 0, 100);
	RooExtendPdf ext_bkg0("ext_bkg0", "extend", bkg0, num_bkg0);
	model0_list.add(ext_bkg0);
	// model0
	RooAddPdf model0("model0", "model", model0_list);


	// 10Be 3.5 MeV
	// RooRealVar x("ex", "ex1", 12.0, 27.0);
	RooDataHist hist1("dh1", "excited energy", x, hist_ex+1);
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
			new RooRealVar("num1_172", "num", 5, 0.1, 20)
		}
	));
	// // 17.7
	// info_map[1].insert(std::make_pair(
	// 	"177",
	// 	VoigtInfo{
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma1_177", "sigma", 0.12*resolution_fraction),
    //         threshold[1],
			// new RooRealVar("num1_177", "num", 5, 0.1, 20)
	// 	}
	// ));
	// 18.5
	info_map[1].insert(std::make_pair(
		"185",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_185", "sigma", 0.13*resolution_fraction),
            threshold[1],
			new RooRealVar("num1_185", "num", 10, 0.1, 50)
		}
	));
	// 19.0
	info_map[1].insert(std::make_pair(
        "190",
        VoigtInfo {
			nullptr,
			nullptr,
            new RooConstVar("sigma1_190", "sigma", 0.15*resolution_fraction),
            threshold[1],
			new RooRealVar("num1_190", "num", 20, 0.1, 50)
        }
    ));
	// 19.3
	info_map[1].insert(std::make_pair(
		"193",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_193", "sigma", 0.16*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_193", "num", 40, 0.1, 100)
		}
	));
	// 20.1
	info_map[1].insert(std::make_pair(
		"201",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_201", "sigma", 0.17*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_201", "num", 30, 0.1, 100)
		}
	));
	// 20.5
	info_map[1].insert(std::make_pair(
		"205",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_205", "sigma", 0.18*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_205", "num", 60, 0.1, 200)
		}
	));
	// 20.9
	info_map[1].insert(std::make_pair(
		"209",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_209", "sigma", 0.18*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_209", "num", 15, 0.1, 100)
		}
	));
	// 216
	info_map[1].insert(std::make_pair(
		"216",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_216", "sigma", 0.19*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_216", "num", 50, 0.1, 200)
		}
	));
	// 223
	info_map[1].insert(std::make_pair(
		"223",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_223", "sigma", 0.20*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_223", "num", 25, 0.1, 100)
		}
	));
	// 228
	info_map[1].insert(std::make_pair(
		"228",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_228", "sigma", 0.21*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_228", "num", 25, 0.1, 100)
		}
	));
	// 232
	info_map[1].insert(std::make_pair(
		"232",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_232", "sigma", 0.21*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_232", "num", 25, 0.1, 100)
		}
	));
	// 236
	info_map[1].insert(std::make_pair(
		"236",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_236", "sigma", 0.21*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_236", "num", 25, 0.1, 100)
		}
	));
	// 245
	info_map[1].insert(std::make_pair(
		"245",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_245", "sigma", 0.22*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_245", "num", 35, 0.1, 100)
		}
	));
	// 254
	info_map[1].insert(std::make_pair(
		"254",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_254", "sigma", 0.23*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_254", "num", 25, 0.1, 100)
		}
	));
	// 260
	info_map[1].insert(std::make_pair(
		"260",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_260", "sigma", 0.25*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_260", "num", 25, 0.1, 100)
		}
	));
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
	// // numebr of events
	// std::vector<RooRealVar*> num_events1;
	// // extend P.D.F.
	// std::vector<RooExtendPdf*> extend_pdf1;
	// argument list for model1
	RooArgList model1_list;
	// construct extend P.D.F.
	for (size_t i = 0; i < pdf1.size(); ++i) {
		// num_events1.push_back(new RooRealVar(
		// 	TString::Format("num1_%ld", i), "num",
		// 	30, 0, 1000
		// ));
		// extend_pdf1.push_back(new RooExtendPdf(
        //     TString::Format("extend1_%ld", i), "extend",
        //     *pdf1[i], *num_events1[i]
        // ));
		model1_list.add(*extend_pdf1[i]);
	}
	// background
	RooRealVar bkg1_p1("bkg1_p1", "p1", 46, 35, 50);
	RooRealVar bkg1_p0("bkg1_p0", "p0", -450, -1000, 0.0);
	BackgroundPoly bkg1("bkg1", "bkg", x, bkg1_p1, bkg1_p0);
	// background events' number
	RooRealVar num_bkg1("num_bkg1", "num", 10, 0, 100);
	RooExtendPdf ext_bkg1("ext_bkg1", "extend", bkg1, num_bkg1);
	model1_list.add(ext_bkg1);
	// model1
	RooAddPdf model1("model1", "model", model1_list);


	// 10Be 6 MeV
	RooDataHist hist2("dh2", "excited energy", RooArgList(x), hist_ex+2);
	// RooDataSet set2("s2", "excited energy", fit_tree[2], RooArgSet(x));
	// construct P.D.F. information
	// 20.1
	info_map[2].insert(std::make_pair(
		"201",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_201", "sigma", 0.11*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_201", "num", 15, 0.1, 100)
		}
	));
	// 20.5
	info_map[2].insert(std::make_pair(
		"205",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_205", "sigma", 0.12*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_205", "num", 50, 0.1, 150)
		}
	));
	// 20.9
	info_map[2].insert(std::make_pair(
		"209",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_209", "sigma", 0.13*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_209", "num", 40, 0.1, 150)
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
			new RooRealVar("num2_216", "num", 50, 0.1, 150)
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
			new RooRealVar("num2_223", "num", 120, 0.1, 300)
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
			new RooRealVar("num2_228", "num", 10, 0.1, 200)
		}
	));
	// 23.2
	info_map[2].insert(std::make_pair(
		"232",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_232", "sigma", 0.17*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_232", "num", 50, 0.1, 200)
		}
	));
	// 23.6
	info_map[2].insert(std::make_pair(
		"236",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_236", "sigma", 0.18*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_236", "num", 50, 0.1, 200)
		}
	));
	// 24.5
	info_map[2].insert(std::make_pair(
		"245",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_245", "sigma", 0.19*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_245", "num", 50, 0.1, 200)
		}
	));
	// 25.4
	info_map[2].insert(std::make_pair(
		"254",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_254", "sigma", 0.20*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_254", "num", 30, 0.1, 200)
		}
	));
	// 26.0
	info_map[2].insert(std::make_pair(
		"260",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_260", "sigma", 0.21*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_260", "num", 30, 0.1, 200)
		}
	));
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
	// // number of events
	// std::vector<RooRealVar*> num_events2;
	// // extended P.D.F.
	// std::vector<RooExtendPdf*> extend_pdf2;
	// arguments list for model0
	RooArgList model2_list;
	// construct extended P.D.F.
	for (size_t i = 0; i < pdf2.size(); ++i) {
		// num_events2.push_back(new RooRealVar(
		// 	TString::Format("num2_%ld", i), "num",
		// 	40, 0, 1000
		// ));
		// extend_pdf2.push_back(new RooExtendPdf(
        //     TString::Format("extend2_%ld", i), "extend",
        //     *pdf2[i], *num_events2[i]
        // ));
		model2_list.add(*extend_pdf2[i]);
	}
	// background
	RooRealVar bkg2_p1("bkg2_p1", "p1", 46, 40, 60);
	RooRealVar bkg2_p0("bkg2_p0", "p0", -520, -1000, 0.0);
	BackgroundPoly bkg2("bkg2", "bkg", x, bkg2_p1, bkg2_p0);
	// background events' number
	RooRealVar num_bkg2("num_bkg2", "num", 10, 0, 100);
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

	// show bar
	std::vector<std::string> bar_names = {
		"172", "185", "190", "193",
		"201", "205", "216", "223",
		"228", "232", "236", "245"
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
	const size_t print_n = 45;
	int print_state[print_n] = {
		0, 0, 0, 0,
		0, 1,
		0, 1,
		0, 1,
		0, 1,
		0, 1, 2,
		0, 1, 2,
		0, 1, 2,
		0, 1, 2,
		0, 1, 2,
		0, 1, 2,
		0, 1, 2,
		0, 1, 2,
		0, 1, 2,
		0, 1, 2,
		0, 1, 2,
	};
	std::string print_name[print_n] = {
		"136", "149", "154", "165",
		"172", "172",
		"185", "185",
		"190", "190",
		"193", "193",
		"201", "201", "201",
		"205", "205", "205",
		"209", "209", "209",
		"216", "216", "216",
		"223", "223", "223",
		"228", "228", "228",
		"232", "232", "232",
		"236", "236", "236",
		"245", "245", "245",
		"254", "254", "254",
		"260", "260", "260"
	};
	for (size_t i = 0; i < print_n; ++i) {
		int state = print_state[i];
		auto search = info_map[state].find(print_name[i]);
		if (search == info_map[state].end()) {
			std::cerr << "Error: Could not find "
				<< print_name[i] << " in state " << state << "\n";
			continue;
		}
		auto &info = search->second;
		double mean = info.mean->getVal();
		double g = info.gamma->getVal();
		double thres = info.threshold;
		double gamma = g * sqrt(mean - thres);
		double event = info.num->getVal();
		double relative_event = event / g_efficiency[state]->Eval(mean);
		std::cout << std::setw(6) << print_state[i]
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
		RooFit::LineColor(kOrange),
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
		RooFit::LineColor(kOrange),
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
		RooFit::LineColor(kOrange),
		RooFit::Components(bkg2)
	);

	// printable efficiency
	TGraph print_efficiency[3];
	const double efficiency_factor[3] = {
		180.0, 300.0, 400.0
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
	// gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	frame0->GetYaxis()->SetLabelSize(0.15);
	frame0->Draw();
	hist_ex[0].GetYaxis()->SetLabelSize(0.15);
	hist_ex[0].SetLineColor(kBlack);
	hist_ex[0].Draw("same");
	print_efficiency[0].Draw("same");
	pads[1]->cd();
	frame1->Draw();
	frame1->GetYaxis()->SetLabelSize(0.15);
	hist_ex[1].SetLineColor(kBlack);
	hist_ex[1].GetYaxis()->SetLabelSize(0.15);
	hist_ex[1].Draw("same");
	print_efficiency[1].Draw("same");
	pads[2]->cd();
	frame2->Draw();
	frame2->GetXaxis()->SetLabelSize(0.12);
	frame2->GetYaxis()->SetLabelSize(0.15);
	hist_ex[2].SetLineColor(kBlack);
	hist_ex[2].GetXaxis()->SetLabelSize(0.12);
	hist_ex[2].GetYaxis()->SetLabelSize(0.15);
	hist_ex[2].Draw("same");
	print_efficiency[2].Draw("same");
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

	TCanvas *c3 = new TCanvas("c3", "", 1920, 1080);
	gStyle->SetOptStat(0);
	TAxis *xaxis = relative_bar[2].GetXaxis();
	for (size_t i = 0; i < bar_names.size(); i++) {
		// get same mean and gamma values
		auto search = info_map[0].find(bar_names[i]);
		if (search == info_map[0].end()) continue;
		double mean = search->second.mean->getVal();
		xaxis->SetNdivisions(13);
		xaxis->ChangeLabel(
			i+2, -1, -1, -1, -1, -1,
			TString::Format("%.1lf", mean)
		);
	}
	relative_bar[2].Draw("bar");
	relative_bar[1].Draw("bar, same");
	relative_bar[0].Draw("bar, same");
	TLegend *legend = new TLegend(0.72, 0.75, 0.88, 0.85);
	legend->AddEntry(relative_bar+0, "ground state", "f");
	legend->AddEntry(relative_bar+1, "2+ state", "f");
	legend->AddEntry(relative_bar+2, "~6MeV state", "f");
	legend->Draw();
	c3->SaveAs(
		TString::Format("%simage/relative-events.png", kGenerateDataPath)
	);
	c3->Write("barc");

	opf.Close();
	ipf.Close();

	// clean

	return 0;
}