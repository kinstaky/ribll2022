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

struct SpectrumV2Event {
	int ppac_flag, taf_flag, target_flag, bind, hole, straight[2];
	double c_kinetic[4], q[4];
	int be_state[4];
	double excited_energy_target[4];
	double stateless_excited_energy[3][4];
};


struct SpectrumV3Event {
	int valid, ppac_flag, taf_flag;
	int be_state;
	double excited_energy;
};


int FillFromSpectrumV2(
	TH1F *hist,
	TH1F *hist_ppac,
	double *tree_ex,
	TTree **tree
) {
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
	SpectrumV2Event event;
	// setup input branches
	ipt->SetBranchAddress("ppac_flag", &event.ppac_flag);
	ipt->SetBranchAddress("taf_flag", &event.taf_flag);
	ipt->SetBranchAddress("target_flag", &event.target_flag);
	ipt->SetBranchAddress("bind", &event.bind);
	ipt->SetBranchAddress("hole", &event.hole);
	ipt->SetBranchAddress("straight", event.straight);
	ipt->SetBranchAddress("c_kinetic", event.c_kinetic);
	ipt->SetBranchAddress("q", event.q);
	ipt->SetBranchAddress("be_state", event.be_state);
	ipt->SetBranchAddress("excited_energy_target", event.excited_energy_target);
	ipt->SetBranchAddress("stateless_excited_energy", event.stateless_excited_energy);

	// loop
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);
		if (event.ppac_flag  == 0) continue;
		if (event.taf_flag != 0) continue;
		if ((event.target_flag & 1) == 0) continue;
		if (event.bind != 0) continue;
		if (event.hole > 0) continue;
		// straight PID
		if (event.straight[0] != 3) continue;
		// beam kinetic energy
		if (event.c_kinetic[0] < 360.0) continue;

		double ex = -1.0;
		int index = -1;
		if (event.q[calc_index] < -11 && event.q[calc_index] > -13) {
			ex = event.stateless_excited_energy[0][calc_index];
			index = 0;
		} else if (event.q[calc_index] < -14.5 && event.q[calc_index] > -16) {
			ex = event.stateless_excited_energy[1][calc_index];
			index = 1;
		} else if (event.q[calc_index] < -17 && event.q[calc_index] > -20) {
			ex = event.stateless_excited_energy[2][calc_index];
			index = 2;
		}
		if (ex > 0) {
			hist[index].Fill(ex);
			tree_ex[index] = ex;
			tree[index]->Fill();
			if (event.ppac_flag == 2) {
				hist_ppac[index].Fill(ex);
			}
		}
	}

	// close file
	ipf.Close();
	return 0;
}


int FillFromSpectrumV3(
	TH1F *hist,
	TH1F *hist_ppac,
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
	ipt->SetBranchAddress("ppac_flag", &event.ppac_flag);
	ipt->SetBranchAddress("taf_flag", &event.taf_flag);
	ipt->SetBranchAddress("be_state", &event.be_state);
	ipt->SetBranchAddress("excited_energy_target", &event.excited_energy);

	// loop
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		if (event.valid != 0) continue;
		if (event.taf_flag != 1) continue;

		if (event.be_state >= 0) {
			hist[event.be_state].Fill(event.excited_energy);
			tree_ex[event.be_state] = event.excited_energy;
			tree[event.be_state]->Fill();
			if (event.ppac_flag == 2) {
				hist_ppac[event.be_state].Fill(event.excited_energy);
			}
		}
	}

	ipf.Close();
	return 0;
}

int main() {	// efficiency file name
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
	// histogram of all valid PPAC events
	TH1F hist_ex[3];
	hist_ex[0] = TH1F("h0", "excited energy", 50, 12.12, 27.12);
	hist_ex[1] = TH1F("h1", "excited energy", 50, 12.12, 27.12);
	hist_ex[2] = TH1F("h2", "excited energy", 50, 12.12, 27.12);
	TH1F hist_ppac_ex[3] = {
		TH1F{"hp0", "excited energy (multiple PPAC)", 50, 12.12, 27.12},
		TH1F{"hp1", "excited energy (multiple PPAC)", 50, 12.12, 27.12},
		TH1F{"hp2", "excited energy (multiple PPAC)", 50, 12.12, 27.12}
	};
	// 分辨率
	TGraph graph_resolution[3];
	// 分支比
	TH1F relative_bar[3];
	relative_bar[0] = TH1F("hb0", "relative events", 15, 0, 15);
	relative_bar[0].SetBarOffset(0.55);
	relative_bar[0].SetBarWidth(0.3);
	relative_bar[0].SetFillColor(kBlue);
	relative_bar[1] = TH1F("hb1", "relative events", 15, 0, 15);
	relative_bar[1].SetBarOffset(0.85);
	relative_bar[1].SetBarWidth(0.3);
	relative_bar[1].SetFillColor(kGreen);
	relative_bar[2] = TH1F("hb2", "relative events", 15, 0, 15);
	relative_bar[2].SetBarOffset(1.15);
	relative_bar[2].SetBarWidth(0.3);
	relative_bar[2].SetFillColor(kRed);

	TTree *fit_tree[3];
	double fit_ex[3];
	for (int i = 0; i < 3; ++i) {
		fit_tree[i] = new TTree("tree", "tree for RooFit");
		fit_tree[i]->Branch("ex", fit_ex+i, "ex/D");
	}

	if (FillFromSpectrumV2(hist_ex, hist_ppac_ex, fit_ex, fit_tree) != 0) {
		return -1;
	}

	opf.cd();
	for (int i = 0; i < 3; ++i) {
		hist_ex[i].Write(TString::Format("ho%d", i));
	}
	for (int i = 0; i < 3; ++i) {
		hist_ppac_ex[i].Write(TString::Format("hp%d", i));
	}

	if (FillFromSpectrumV3(hist_ex, hist_ppac_ex, fit_ex, fit_tree) != 0) {
		return -1;
	}
	opf.cd();
	for (int i = 0; i < 3; ++i) {
		hist_ex[i].Write(TString::Format("hoe%d", i));
	}
	for (int i = 0; i < 3; ++i) {
		hist_ppac_ex[i].Write(TString::Format("hpe%d", i));
	}

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
	RooRealVar x("ex", "excited energy", 12.12, 27.12);

	RooDataHist hist0("dh0", "excited energy", x, hist_ex);
	// RooDataSet set0("s0", "excited energy", fit_tree[0], RooArgSet(x));
	// construct P.D.F. information
	// 13.6
	// info_map[0].insert(std::make_pair(
	// 	"136",
	// 	VoigtInfo {
	// 		new RooRealVar("mean136", "mean", 13.65, 13.6, 13.7),
	// 		new RooRealVar("gamma136", "gamma", 4e-6, 1e-6, 5e-6),
	// 		new RooConstVar("sigma0_136", "sigma", 0.1*resolution_fraction),
	// 		threshold[0],
	// 		new RooRealVar("num0_136", "num", 3, 1, 5)
	// 	}
	// ));
	// 14.4
	info_map[0].insert(std::make_pair(
        "144",
        VoigtInfo {
            new RooRealVar("mean144", "mean", 14.4, 14.3, 14.5),
            new RooRealVar("gamma144", "gamma", 3e-6, 1e-6, 5e-6),
            new RooConstVar("sigma0_144", "sigma", 0.12*resolution_fraction),
            threshold[0],
			new RooRealVar("num0_144", "num", 9, 8, 10)
        }
    ));
	// 14.9
	info_map[0].insert(std::make_pair(
        "149",
        VoigtInfo {
            new RooRealVar("mean149", "mean", 14.9, 14.8, 15.0),
            new RooRealVar("gamma149", "gamma", 0.059, 0.02, 0.07),
            new RooConstVar("sigma0_149", "sigma", 0.14*resolution_fraction),
            threshold[0],
			new RooRealVar("num0_149", "num", 32, 10, 35)
        }
    ));
	// 15.6
	info_map[0].insert(std::make_pair(
        "156",
        VoigtInfo {
            new RooRealVar("mean154", "mean", 15.5, 15.6, 15.7),
            new RooRealVar("gamma154", "gamma", 0.095, 0.07, 0.12),
            new RooConstVar("sigma0_154", "sigma", 0.14*resolution_fraction),
            threshold[0],
			new RooRealVar("num0_154", "num", 83, 50, 85)
        }
    ));
    // 16.5
    info_map[0].insert(std::make_pair(
        "165",
        VoigtInfo {
            new RooRealVar("mean165", "mean", 16.3, 16.4, 16.5),
			new RooRealVar("gamma165", "gamma", 0.06, 0.05, 0.07),
            new RooConstVar("sigma0_165", "sigma", 0.16*resolution_fraction),
            threshold[0],
			new RooRealVar("num0_165", "num", 80, 50, 90)
		}
	));
	// 17.2
	info_map[0].insert(std::make_pair(
        "172",
        VoigtInfo {
            new RooRealVar("mean172", "mean", 17.2, 17.3, 17.4),
            new RooRealVar("gamma172", "gamma", 0.052, 0.04, 0.065),
            new RooConstVar("sigma0_172", "sigma", 0.18*resolution_fraction),
            threshold[1],
			new RooRealVar("num0_172", "num", 10, 0.01, 150)
        }
    ));
    // // 17.5
    // info_map[0].insert(std::make_pair(
    //     "175",
    //     VoigtInfo {
    //         new RooRealVar("mean175", "mean", 17.5, 17.4, 17.8),
    //         new RooRealVar("gamma175", "gamma", 0.2, 1e-6, 0.4),
	// 		new RooConstVar("sigma0_175", "sigma", 0.19*resolution_fraction),
    //         threshold[1],
	// 		new RooRealVar("num0_175", "num", 50, 30, 80)
	// 	}
	// ));
	// 18.3
	info_map[0].insert(std::make_pair(
        "183",
        VoigtInfo {
            new RooRealVar("mean183", "mean", 18.25, 18.2, 18.6),
            new RooRealVar("gamma183", "gamma", 0.05, 1e-6, 0.2),
            new RooConstVar("sigma0_183", "sigma", 0.2*resolution_fraction),
            threshold[1],
			new RooRealVar("num0_183", "num", 115, 40, 150)
        }
    ));
	// 19.0
	info_map[0].insert(std::make_pair(
        "192",
        VoigtInfo {
            new RooRealVar("mean192", "mean", 19.2, 19.0, 19.4),
            new RooRealVar("gamma192", "gamma", 0.2, 0.05, 0.5),
            new RooConstVar("sigma0_192", "sigma", 0.2*resolution_fraction),
            threshold[1],
			new RooRealVar("num0_192", "num", 100, 1.0, 200)
        }
    ));
	// 19.8
	info_map[0].insert(std::make_pair(
        "198",
        VoigtInfo {
            new RooRealVar("mean198", "mean", 19.7, 19.65, 20.1),
            new RooRealVar("gamma198", "gamma", 1e-4, 1e-6, 0.1),
            new RooConstVar("sigma0_198", "sigma", 0.21*resolution_fraction),
            threshold[1],
			new RooRealVar("num0_198", "num", 10, 2, 20)
        }
    ));
	// 20.2
	info_map[0].insert(std::make_pair(
        "202",
        VoigtInfo {
            new RooRealVar("mean202", "mean", 20.4, 20.2, 20.4),
            new RooRealVar("gamma202", "gamma", 0.04, 1e-6, 0.1),
            new RooConstVar("sigma0_202", "sigma", 0.22*resolution_fraction),
            threshold[2],
			new RooRealVar("num0_202", "num", 50, 20, 100)
        }
    ));
    // 20.7
    info_map[0].insert(std::make_pair(
        "207",
        VoigtInfo {
            new RooRealVar("mean207", "mean", 20.6, 20.5, 20.7),
            new RooRealVar("gamma207", "gamma", 0.2, 0.12, 0.3),
			new RooConstVar("sigma0_207", "sigma", 0.23*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_207", "num", 1, 0.01, 30)
		}
	));
	// 21.6
	info_map[0].insert(std::make_pair(
		"216",
		VoigtInfo {
            new RooRealVar("mean216", "mean", 21.4, 21.3, 21.7),
            new RooRealVar("gamma216", "gamma", 0.12, 0.1, 0.5),
            new RooConstVar("sigma0_216", "sigma", 0.24*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_216", "num", 50, 0.01, 100)
		}
	));
	// 22.3
	info_map[0].insert(std::make_pair(
		"223",
		VoigtInfo {
            new RooRealVar("mean223", "mean", 22.25, 22.2, 22.3),
            new RooRealVar("gamma223", "gamma", 0.6, 0.2, 0.8),
            new RooConstVar("sigma0_223", "sigma", 0.25*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_223", "num", 10, 0.01, 40)
		}
	));
	// 22.8
	info_map[0].insert(std::make_pair(
		"228",
		VoigtInfo {
            new RooRealVar("mean228", "mean", 22.7, 22.7, 22.8),
            new RooRealVar("gamma228", "gamma", 0.05, 1e-6, 0.1),
            new RooConstVar("sigma0_228", "sigma", 0.26*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_228", "num", 20, 0.01, 50)
		}
	));
	// 23.2
	info_map[0].insert(std::make_pair(
		"232",
		VoigtInfo {
            new RooRealVar("mean232", "mean", 23.1, 23.0, 23.3),
            new RooRealVar("gamma232", "gamma", 1e-5, 1e-6, 1e-3),
            new RooConstVar("sigma0_232", "sigma", 0.26*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_232", "num", 10, 0.01, 50)
		}
	));
	// 23.7
	info_map[0].insert(std::make_pair(
		"237",
		VoigtInfo {
            new RooRealVar("mean237", "mean", 23.7, 23.6, 23.8),
            new RooRealVar("gamma237", "gamma", 0.1, 0.01, 0.5),
            new RooConstVar("sigma0_237", "sigma", 0.27*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_237", "num", 20, 0.01, 50)
		}
	));
	// 24.0
	info_map[0].insert(std::make_pair(
		"240",
		VoigtInfo {
            new RooRealVar("mean240", "mean", 24.0, 23.9, 24.1),
            new RooRealVar("gamma240", "gamma", 1e-6, 1e-6, 1e-3),
            new RooConstVar("sigma0_240", "sigma", 0.27*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_240", "num", 5, 0.01, 20)
		}
	));
	// 24.7
	info_map[0].insert(std::make_pair(
		"247",
		VoigtInfo {
            new RooRealVar("mean247", "mean", 24.5, 24.3, 24.8),
            new RooRealVar("gamma247", "gamma", 0.01, 1e-6, 0.2),
            new RooConstVar("sigma0_247", "sigma", 0.29*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_247", "num", 1, 0.01, 10)
		}
	));
	// // 25.4
	// info_map[0].insert(std::make_pair(
	// 	"254",
	// 	VoigtInfo {
    //         new RooRealVar("mean254", "mean", 25.4, 25.2, 25.7),
    //         new RooRealVar("gamma254", "gamma", 0.02, 0.0, 1.0),
    //         new RooConstVar("sigma0_254", "sigma", 0.28*resolution_fraction),
	// 		threshold[2],
	// 		new RooRealVar("num0_254", "num", 1, 0.01, 10)
	// 	}
	// ));
	// 26.0
	info_map[0].insert(std::make_pair(
		"260",
		VoigtInfo {
            new RooRealVar("mean260", "mean", 26.1, 26.0, 26.2),
            new RooRealVar("gamma260", "gamma", 0.5, 1e-6, 2.0),
            new RooConstVar("sigma0_260", "sigma", 0.31*resolution_fraction),
			threshold[2],
			new RooRealVar("num0_260", "num", 4, 0.01, 20)
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
	RooRealVar bkg0_root0("bkg0_root0", "root0", 14.0, 14.0, 16.0);
	RooRealVar bkg0_root1("bkg0_root1", "root1", 26.0, 25.0, 27.0);
	BackgroundPoly bkg0("bkg0", "bkg", x, bkg0_root0, bkg0_root1);
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
			new RooRealVar("num1_172", "num", 30, 20, 150)
		}
	));
	// 17.5
	// info_map[1].insert(std::make_pair(
	// 	"175",
	// 	VoigtInfo{
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma1_175", "sigma", 0.12*resolution_fraction),
    //         threshold[1],
	// 		new RooRealVar("num1_175", "num", 5, 0.01, 50)
	// 	}
	// ));
	// 18.5
	info_map[1].insert(std::make_pair(
		"183",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_183", "sigma", 0.13*resolution_fraction),
            threshold[1],
			new RooRealVar("num1_183", "num", 100, 20, 150)
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
	// 20.2
	info_map[1].insert(std::make_pair(
		"202",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_202", "sigma", 0.17*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_202", "num", 20, 2, 30)
		}
	));
	// // 20.5
	// info_map[1].insert(std::make_pair(
	// 	"205",
	// 	VoigtInfo{
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma1_205", "sigma", 0.18*resolution_fraction),
    //         threshold[2],
	// 		new RooRealVar("num1_205", "num", 60, 0.1, 200)
	// 	}
	// ));
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
			new RooRealVar("num1_223", "num", 5, 0.1, 100)
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
			new RooRealVar("num1_228", "num", 70, 5, 140)
		}
	));
	// 23.2
	info_map[1].insert(std::make_pair(
		"232",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_232", "sigma", 0.21*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_232", "num", 10, 0.1, 50)
		}
	));
	// 23.7
	info_map[1].insert(std::make_pair(
		"237",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_237", "sigma", 0.21*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_237", "num", 30, 10, 70)
		}
	));
	// 24.0
	info_map[1].insert(std::make_pair(
		"240",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_240", "sigma", 0.23*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_240", "num", 25, 0.1, 100)
		}
	));
	// 24.7
	info_map[1].insert(std::make_pair(
		"247",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_247", "sigma", 0.23*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_247", "num", 25, 0.1, 100)
		}
	));
	// 254
	// info_map[1].insert(std::make_pair(
	// 	"254",
	// 	VoigtInfo{
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma1_254", "sigma", 0.23*resolution_fraction),
    //         threshold[2],
	// 		new RooRealVar("num1_254", "num", 25, 0.1, 100)
	// 	}
	// ));
	// 260
	info_map[1].insert(std::make_pair(
		"260",
		VoigtInfo{
			nullptr,
			nullptr,
			new RooConstVar("sigma1_260", "sigma", 0.25*resolution_fraction),
            threshold[2],
			new RooRealVar("num1_260", "num", 100, 0.1, 200)
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
	RooRealVar bkg1_root0("bkg1_root0", "root0", 17, 16, 18);
	RooRealVar bkg1_root1("bkg1_root1", "root1", 26, 25, 27);
	BackgroundPoly bkg1("bkg1", "bkg", x, bkg1_root0, bkg1_root1);
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
	// 20.2
	info_map[2].insert(std::make_pair(
		"202",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_202", "sigma", 0.11*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_202", "num", 5, 1, 20)
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
			new RooRealVar("num2_228", "num", 10, 0.1, 50)
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
			new RooRealVar("num2_232", "num", 50, 20, 100)
		}
	));
	// 23.6
	info_map[2].insert(std::make_pair(
		"237",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_237", "sigma", 0.18*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_237", "num", 1, 0.01, 10)
		}
	));
	// 24.0
	info_map[2].insert(std::make_pair(
		"240",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_240", "sigma", 0.18*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_240", "num", 60, 40, 100)
		}
	));
	// 24.7
	info_map[2].insert(std::make_pair(
		"247",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_247", "sigma", 0.20*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_247", "num", 50, 0.1, 200)
		}
	));
	// // 25.4
	// info_map[2].insert(std::make_pair(
	// 	"254",
	// 	VoigtInfo {
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma2_254", "sigma", 0.20*resolution_fraction),
	// 		threshold[2],
	// 		new RooRealVar("num2_254", "num", 30, 0.1, 200)
	// 	}
	// ));
	// 26.0
	info_map[2].insert(std::make_pair(
		"260",
		VoigtInfo {
			nullptr,
			nullptr,
			new RooConstVar("sigma2_260", "sigma", 0.21*resolution_fraction),
			threshold[2],
			new RooRealVar("num2_260", "num", 100, 0.1, 200)
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
	RooRealVar bkg2_root0("bkg2_root0", "root0", 20.0, 19.5, 21.0);
	RooRealVar bkg2_root1("bkg2_root1", "root1", 26.0, 25.0, 27.0);
	BackgroundPoly bkg2("bkg2", "bkg", x, bkg2_root0, bkg2_root1);
	// background events' number
	RooRealVar num_bkg2("num_bkg2", "num", 10, 0, 50);
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
		"172", "175", "183", "192",
		"198", "202", "207", "216",
		"223", "228", "232", "237",
		"240", "247"
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
		{0, "136"},
		{0, "144"},
		{0, "149"},
		{0, "156"},
		{0, "165"},
		{0, "172"},
		{1, "172"},
		{0, "175"},
		{1, "175"},
		{0, "183"},
		{1, "183"},
		{0, "192"},
		{1, "192"},
		{0, "198"},
		{1, "198"},
		{2, "198"},
		{0, "202"},
		{1, "202"},
		{2, "202"},
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
		{0, "232"},
		{1, "232"},
		{2, "232"},
		{0, "237"},
		{1, "237"},
		{2, "237"},
		{0, "240"},
		{1, "240"},
		{2, "240"},
		{0, "247"},
		{1, "247"},
		{2, "247"},
		{0, "260"},
		{1, "260"},
		{2, "260"},
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
		250.0, 450.0, 400.0
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
	for (const auto &[key, value] : info_map[0]) {
		TLine *line = new TLine(
			value.mean->getVal(), 0.0, value.mean->getVal(), 35.0
		);
		line->SetLineStyle(2);
		line->Draw("same");
	}
	print_efficiency[0].Draw("same");

	pads[1]->cd();
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
		TString::Format("%s%sroofit-result.png", kGenerateDataPath, kImageDir)
	);
	c1->Write("c1");

	TCanvas *c2 = new TCanvas("c2", "", 1920, 1080);
	test_frame->Draw();
	c2->SaveAs(
		TString::Format("%s%stest-asym-voigt.png", kGenerateDataPath, kImageDir)
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
		xaxis->SetNdivisions(bar_names.size()+1);
		xaxis->ChangeLabel(
			i+2, -1, -1, -1, -1, -1,
			TString::Format("%.1lf", mean)
		);
	}
	xaxis->ChangeLabel(1, -1, 0, -1, -1, -1, "");
	xaxis->ChangeLabel(bar_names.size()+2, -1, 0, -1, -1, -1, "");
	relative_bar[2].Draw("bar");
	relative_bar[1].Draw("bar, same");
	relative_bar[0].Draw("bar, same");
	TLegend *legend = new TLegend(0.72, 0.75, 0.88, 0.85);
	legend->AddEntry(relative_bar+0, "ground state", "f");
	legend->AddEntry(relative_bar+1, "2+ state", "f");
	legend->AddEntry(relative_bar+2, "~6MeV state", "f");
	legend->Draw();
	c3->SaveAs(
		TString::Format("%s%srelative-events.png", kGenerateDataPath, kImageDir)
	);
	c3->Write("barc");

	opf.Close();

	return 0;
}