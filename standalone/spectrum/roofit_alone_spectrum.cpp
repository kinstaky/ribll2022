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
constexpr double resolution_fraction[3] = {1.4, 1.4, 1.4};

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
	const RooAbsPdf &component
) {
	pdf.plotOn(
		frame,
		RooFit::LineColor(kBlue),
		RooFit::Components(component)
	);
}

struct SpectrumV2Event {
	int ppac_flag, taf_flag, target_flag, bind, hole, straight[2];
	double q[4];
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
	TH1F *hist
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

		double ex = -1.0;
		int index = -1;
		if (event.q[0] < -11 && event.q[0] > -13) {
			ex = event.stateless_excited_energy[0][0];
			index = 0;
		} else if (event.q[0] < -14.5 && event.q[0] > -16) {
			ex = event.stateless_excited_energy[1][0];
			index = 1;
		} else if (event.q[0] < -17 && event.q[0] > -20) {
			ex = event.stateless_excited_energy[2][0];
			index = 2;
		}
		if (ex > 0) {
			hist[index].Fill(ex);
		}
	}

	// close file
	ipf.Close();
	return 0;
}


int FillFromSpectrumV3(
	TH1F *hist
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
	// 分辨率
	TGraph graph_resolution[3];
	// 分支比
	// TH1F relative_bar[3];
	// relative_bar[0] = TH1F("hb0", "relative events", 13, 0, 13);
	// relative_bar[0].SetBarOffset(0.55);
	// relative_bar[0].SetBarWidth(0.3);
	// relative_bar[0].SetFillColor(kBlue);
	// relative_bar[1] = TH1F("hb1", "relative events", 13, 0, 13);
	// relative_bar[1].SetBarOffset(0.85);
	// relative_bar[1].SetBarWidth(0.3);
	// relative_bar[1].SetFillColor(kGreen);
	// relative_bar[2] = TH1F("hb2", "relative events", 13, 0, 13);
	// relative_bar[2].SetBarOffset(1.15);
	// relative_bar[2].SetBarWidth(0.3);
	// relative_bar[2].SetFillColor(kRed);


	// fill events
	if (FillFromSpectrumV2(hist_ex) != 0) {
		return -1;
	}
	if (FillFromSpectrumV3(hist_ex) != 0) {
		return -1;
	}


	// fitting function information
	std::map<std::string, VoigtInfo> info_map[3];

	// RooFit
	// get data
	RooRealVar x0("ex0", "excited energy 0", 12.12, 27.12);
	RooRealVar x1("ex1", "excited energy 1", 12.12, 27.12);
	RooRealVar x2("ex2", "excited energy 2", 12.12, 27.12);

	RooDataHist hist0("dh0", "excited energy", x0, hist_ex);
	// construct P.D.F. information
	// 13.6
	info_map[0].insert(std::make_pair(
		"136",
		VoigtInfo {
			new RooRealVar("mean0_136", "mean", 13.65, 13.6, 13.7),
			new RooRealVar("gamma0_136", "gamma", 0.05, 0.0, 0.2),
			new RooConstVar("sigma0_136", "sigma", 0.1*resolution_fraction[0]),
			threshold[0],
			new RooRealVar("num0_136", "num", 3, 0.01, 20)
		}
	));
	// 14.4
	info_map[0].insert(std::make_pair(
        "144",
        VoigtInfo {
            new RooRealVar("mean0_144", "mean", 14.4, 14.2, 14.7),
            new RooRealVar("gamma0_144", "gamma", 0.06, 0.0, 0.5),
            new RooConstVar("sigma0_144", "sigma", 0.12*resolution_fraction[0]),
            threshold[0],
			new RooRealVar("num0_144", "num", 10, 0.01, 50)
        }
    ));
	// 14.9
	info_map[0].insert(std::make_pair(
        "149",
        VoigtInfo {
            new RooRealVar("mean0_149", "mean", 14.9, 14.8, 15.3),
            new RooRealVar("gamma0_149", "gamma", 0.1, 0.0, 0.4),
            new RooConstVar("sigma0_149", "sigma", 0.14*resolution_fraction[0]),
            threshold[0],
			new RooRealVar("num0_149", "num", 40, 0.01, 100)
        }
    ));
	// 15.6
	info_map[0].insert(std::make_pair(
        "156",
        VoigtInfo {
            new RooRealVar("mean0_156", "mean", 15.6, 15.4, 16.0),
            new RooRealVar("gamma0_156", "gamma", 0.1, 0.0, 0.4),
            new RooConstVar("sigma0_156", "sigma", 0.14*resolution_fraction[0]),
            threshold[0],
			new RooRealVar("num0_156", "num", 40, 0.01, 100)
        }
    ));
    // 16.5
    info_map[0].insert(std::make_pair(
        "165",
        VoigtInfo {
            new RooRealVar("mean0_165", "mean", 16.4, 16.2, 17.0),
			new RooRealVar("gamma0_165", "gamma", 0.1, 0.0, 0.4),
            new RooConstVar("sigma0_165", "sigma", 0.16*resolution_fraction[0]),
            threshold[0],
			new RooRealVar("num0_165", "num", 40, 0.01, 100)
		}
	));
    // 17.2
    info_map[0].insert(std::make_pair(
        "172",
        VoigtInfo {
            new RooRealVar("mean0_172", "mean", 17.2, 17.19, 17.21),
            new RooRealVar("gamma0_172", "gamma", 0.3, 0.29, 0.31),
			new RooConstVar("sigma0_172", "sigma", 0.18*resolution_fraction[0]),
            threshold[1],
			new RooRealVar("num0_172", "num", 20, 0.01, 50)
		}
	));
	// 17.5
	info_map[0].insert(std::make_pair(
        "175",
        VoigtInfo {
            new RooRealVar("mean0_175", "mean", 17.3, 16.9, 17.6),
            new RooRealVar("gamma0_175", "gamma", 0.1, 0.0, 0.4),
            new RooConstVar("sigma0_175", "sigma", 0.18*resolution_fraction[0]),
            threshold[1],
			new RooRealVar("num0_175", "num", 50, 0.01, 100)
        }
    ));
	// 18.3
	info_map[0].insert(std::make_pair(
        "183",
        VoigtInfo {
            new RooRealVar("mean0_183", "mean", 18.3, 18.0, 18.6),
            new RooRealVar("gamma0_183", "gamma", 0.05, 0.0, 0.4),
            new RooConstVar("sigma0_183", "sigma", 0.2*resolution_fraction[0]),
            threshold[1],
			new RooRealVar("num0_183", "num", 80, 0.01, 150)
        }
    ));
	// // 19.0
	// info_map[0].insert(std::make_pair(
    //     "190",
    //     VoigtInfo {
    //         new RooRealVar("mean190", "mean", 19.0, 18.5, 19.2),
    //         new RooRealVar("gamma190", "gamma", 0.05, 0.0, 0.4),
    //         new RooConstVar("sigma0_190", "sigma", 0.2*resolution_fraction[0]),
    //         threshold[1],
	// 		new RooRealVar("num0_190", "num", 1, 0.01, 10)
    //     }
    // ));
    // 19.5
    info_map[0].insert(std::make_pair(
        "195",
        VoigtInfo {
            new RooRealVar("mean0_195", "mean", 19.5, 19.2, 19.7),
            new RooRealVar("gamma0_195", "gamma", 0.05, 0.0, 1.0),
			new RooConstVar("sigma0_195", "sigma", 0.22*resolution_fraction[0]),
			threshold[2],
			new RooRealVar("num0_195", "num", 50, 0.01, 150)
		}
	));
	// 20.2
	info_map[0].insert(std::make_pair(
        "202",
        VoigtInfo {
            new RooRealVar("mean0_202", "mean", 20.2, 19.8, 20.5),
            new RooRealVar("gamma0_202", "gamma", 0.02, 0.0, 0.4),
            new RooConstVar("sigma0_202", "sigma", 0.22*resolution_fraction[0]),
            threshold[2],
			new RooRealVar("num0_201", "num", 20, 0.01, 100)
        }
    ));
    // // 20.5
    // info_map[0].insert(std::make_pair(
    //     "205",
    //     VoigtInfo {
    //         new RooRealVar("mean205", "mean", 20.5, 20.3, 20.7),
    //         new RooRealVar("gamma205", "gamma", 0.02, 0.0, 0.4),
	// 		new RooConstVar("sigma0_205", "sigma", 0.23*resolution_fraction[0]),
	// 		threshold[2],
	// 		new RooRealVar("num0_205", "num", 30, 0.01, 50)
	// 	}
	// ));
	// 20.9
	// info_map[0].insert(std::make_pair(
	// 	"209",
	// 	VoigtInfo {
	// 		new RooRealVar("mean209", "mean", 20.9, 20.7, 21.3),
    //         new RooRealVar("gamma209", "gamma", 0.05, 0.0, 0.2),
    //         new RooConstVar("sigma0_209", "sigma", 0.23*resolution_fraction[0]),
    //         threshold[2],
	// 		new RooRealVar("num0_209", "num", 30, 0.01, 100)
	// 	}
	// ));
	// 21.6
	info_map[0].insert(std::make_pair(
		"216",
		VoigtInfo {
            new RooRealVar("mean0_216", "mean", 21.6, 21.2, 21.8),
            new RooRealVar("gamma0_216", "gamma", 0.05, 0.0, 0.2),
            new RooConstVar("sigma0_216", "sigma", 0.24*resolution_fraction[0]),
			threshold[2],
			new RooRealVar("num0_216", "num", 30, 0.01, 50)
		}
	));
	// 22.3
	info_map[0].insert(std::make_pair(
		"223",
		VoigtInfo {
            new RooRealVar("mean0_223", "mean", 22.3, 22.0, 22.6),
            new RooRealVar("gamma0_223", "gamma", 0.05, 0.0, 1.0),
            new RooConstVar("sigma0_223", "sigma", 0.25*resolution_fraction[0]),
			threshold[2],
			new RooRealVar("num0_223", "num", 5, 0.01, 20)
		}
	));
	// 22.8
	info_map[0].insert(std::make_pair(
		"228",
		VoigtInfo {
            new RooRealVar("mean0_228", "mean", 22.8, 22.7, 23.1),
            new RooRealVar("gamma0_228", "gamma", 0.05, 0.0, 0.2),
            new RooConstVar("sigma0_228", "sigma", 0.26*resolution_fraction[0]),
			threshold[2],
			new RooRealVar("num0_228", "num", 20, 0.01, 50)
		}
	));
	// 23.2
	info_map[0].insert(std::make_pair(
		"232",
		VoigtInfo {
            new RooRealVar("mean0_232", "mean", 23.2, 23.0, 23.5),
            new RooRealVar("gamma0_232", "gamma", 0.02, 0.0, 1.0),
            new RooConstVar("sigma0_232", "sigma", 0.26*resolution_fraction[0]),
			threshold[2],
			new RooRealVar("num0_232", "num", 5, 0.01, 50)
		}
	));
	// 23.6
	info_map[0].insert(std::make_pair(
		"236",
		VoigtInfo {
            new RooRealVar("mean0_236", "mean", 23.6, 23.3, 24.0),
            new RooRealVar("gamma0_236", "gamma", 0.02, 0.0, 0.4),
            new RooConstVar("sigma0_236", "sigma", 0.27*resolution_fraction[0]),
			threshold[2],
			new RooRealVar("num0_236", "num", 10, 0.01, 50)
		}
	));
	// 24.5
	info_map[0].insert(std::make_pair(
		"245",
		VoigtInfo {
            new RooRealVar("mean0_245", "mean", 24.5, 24.2, 24.7),
            new RooRealVar("gamma0_245", "gamma", 0.02, 0.0, 0.8),
            new RooConstVar("sigma0_245", "sigma", 0.28*resolution_fraction[0]),
			threshold[2],
			new RooRealVar("num0_245", "num", 1, 0.01, 10)
		}
	));
	// 25.4
	info_map[0].insert(std::make_pair(
		"254",
		VoigtInfo {
            new RooRealVar("mean0_254", "mean", 25.4, 25.2, 25.7),
            new RooRealVar("gamma0_254", "gamma", 0.02, 0.0, 2.0),
            new RooConstVar("sigma0_254", "sigma", 0.28*resolution_fraction[0]),
			threshold[2],
			new RooRealVar("num0_254", "num", 1, 0.01, 10)
		}
	));
	// 26.0
	info_map[0].insert(std::make_pair(
		"260",
		VoigtInfo {
            new RooRealVar("mean0_260", "mean", 26.2, 26.0, 26.4),
            new RooRealVar("gamma0_260", "gamma", 0.02, 0.0, 1.0),
            new RooConstVar("sigma0_260", "sigma", 0.31*resolution_fraction[0]),
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
			x0,
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
	// // background
	RooRealVar bkg0_root0("bkg0_root0", "root0", 14.0, 14.0, 16.0);
	RooRealVar bkg0_root1("bkg0_root1", "root1", 26.0, 25.0, 27.0);
	BackgroundPoly bkg0("bkg0", "bkg", x0, bkg0_root0, bkg0_root1);
	// background events' number
	RooRealVar num_bkg0("num_bkg0", "num", 10, 0, 300);
	RooExtendPdf ext_bkg0("ext_bkg0", "extend", bkg0, num_bkg0);
	model0_list.add(ext_bkg0);
	// model0
	RooAddPdf model0("model0", "model", model0_list);


	// 10Be 3.5 MeV
	// RooRealVar x("ex", "ex1", 12.0, 27.0);
	RooDataHist hist1("dh1", "excited energy", x1, hist_ex+1);
	// RooDataSet set1("s1", "excited energy", fit_tree[1], RooArgSet(x1));
	// construct P.D.F. information
	// 17.2
	info_map[1].insert(std::make_pair(
		"172",
		VoigtInfo{
			new RooRealVar("mean1_172", "mean", 17.3, 16.9, 17.6),
			new RooRealVar("gamma1_172", "gamma", 0.1, 0.0, 0.4),
			new RooConstVar("sigma1_172", "sigma", 0.11*resolution_fraction[1]),
            threshold[1],
			new RooRealVar("num1_172", "num", 10, 0.1, 200)
		}
	));
	// // 17.7
	// info_map[1].insert(std::make_pair(
	// 	"177",
	// 	VoigtInfo{
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma1_177", "sigma", 0.12*resolution_fraction[1]),
    //         threshold[1],
			// new RooRealVar("num1_177", "num", 5, 0.1, 20)
	// 	}
	// ));
	// 18.3
	info_map[1].insert(std::make_pair(
		"183",
		VoigtInfo{
			new RooRealVar("mean1_183", "mean", 18.3, 18.0, 18.6),
			new RooRealVar("gamma1_183", "gamma", 0.05, 0.0, 1.0),
			new RooConstVar("sigma1_183", "sigma", 0.13*resolution_fraction[1]),
            threshold[1],
			new RooRealVar("num1_183", "num", 50, 0.1, 500)
		}
	));
	// 19.0
	// info_map[1].insert(std::make_pair(
    //     "190",
    //     VoigtInfo {
	// 		nullptr,
	// 		nullptr,
    //         new RooConstVar("sigma1_190", "sigma", 0.15*resolution_fraction[1]),
    //         threshold[1],
	// 		new RooRealVar("num1_190", "num", 20, 0.1, 50)
    //     }
    // ));
	// 19.5
	info_map[1].insert(std::make_pair(
		"195",
		VoigtInfo{
			new RooRealVar("mean1_195", "mean", 19.5, 19.2, 19.7),
			new RooRealVar("gamma1_195", "gamma", 0.05, 0.0, 1.0),
			new RooConstVar("sigma1_195", "sigma", 0.16*resolution_fraction[1]),
            threshold[2],
			new RooRealVar("num1_195", "num", 100, 0.1, 500)
		}
	));
	// // 20.1
	// info_map[1].insert(std::make_pair(
	// 	"201",
	// 	VoigtInfo{
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma1_201", "sigma", 0.17*resolution_fraction[1]),
    //         threshold[2],
	// 		new RooRealVar("num1_201", "num", 30, 0.1, 100)
	// 	}
	// ));
	// // 20.5
	// info_map[1].insert(std::make_pair(
	// 	"205",
	// 	VoigtInfo{
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma1_205", "sigma", 0.18*resolution_fraction[1]),
    //         threshold[2],
	// 		new RooRealVar("num1_205", "num", 60, 0.1, 200)
	// 	}
	// ));
	// 20.7
	info_map[1].insert(std::make_pair(
		"207",
		VoigtInfo{
			new RooRealVar("mean1_207", "mean", 20.7, 20.5, 20.8),
			new RooRealVar("gamma1_207", "gamma", 0.02, 0.0, 1.0),
			new RooConstVar("sigma1_209", "sigma", 0.18*resolution_fraction[1]),
            threshold[2],
			new RooRealVar("num1_209", "num", 50, 0.1, 500)
		}
	));
	// 216
	info_map[1].insert(std::make_pair(
		"216",
		VoigtInfo{
			new RooRealVar("mean1_216", "mean", 21.6, 21.3, 21.8),
			new RooRealVar("gamma1_216", "gamma", 0.05, 0.0, 0.2),
			new RooConstVar("sigma1_216", "sigma", 0.19*resolution_fraction[1]),
            threshold[2],
			new RooRealVar("num1_216", "num", 50, 0.1, 200)
		}
	));
	// 223
	// info_map[1].insert(std::make_pair(
	// 	"223",
	// 	VoigtInfo{
	// 		new RooRealVar("mean223", "mean", 22.3, 22.0, 22.6),
	// 		new RooRealVar("gamma223", "gamma", 0.05, 0.0, 0.2),
	// 		new RooConstVar("sigma1_223", "sigma", 0.20*resolution_fraction[1]),
    //         threshold[2],
	// 		new RooRealVar("num1_223", "num", 25, 0.1, 100)
	// 	}
	// ));
	// 228
	info_map[1].insert(std::make_pair(
		"228",
		VoigtInfo{
			new RooRealVar("mean1_228", "mean", 22.8, 22.7, 23.1),
			new RooRealVar("gamma1_228", "gamma", 0.05, 0.0, 0.2),
			new RooConstVar("sigma1_228", "sigma", 0.21*resolution_fraction[1]),
            threshold[2],
			new RooRealVar("num1_228", "num", 25, 0.1, 100)
		}
	));
	// 232
	// info_map[1].insert(std::make_pair(
	// 	"232",
	// 	VoigtInfo{
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma1_232", "sigma", 0.21*resolution_fraction[1]),
    //         threshold[2],
	// 		new RooRealVar("num1_232", "num", 25, 0.1, 100)
	// 	}
	// ));
	// 236
	info_map[1].insert(std::make_pair(
		"236",
		VoigtInfo{
			new RooRealVar("mean1_236", "mean", 23.6, 23.2, 24.0),
			new RooRealVar("gamma1_236", "gamma", 0.05, 0.0, 0.2),
			new RooConstVar("sigma1_236", "sigma", 0.21*resolution_fraction[1]),
            threshold[2],
			new RooRealVar("num1_236", "num", 25, 0.1, 100)
		}
	));
	// 245
	info_map[1].insert(std::make_pair(
		"245",
		VoigtInfo{
			new RooRealVar("mean1_245", "mean", 24.5, 24.2, 24.7),
			new RooRealVar("gamma1_245", "gamma", 0.05, 0.0, 0.5),
			new RooConstVar("sigma1_245", "sigma", 0.22*resolution_fraction[1]),
            threshold[2],
			new RooRealVar("num1_245", "num", 35, 0.1, 100)
		}
	));
	// 254
	// info_map[1].insert(std::make_pair(
	// 	"254",
	// 	VoigtInfo{
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma1_254", "sigma", 0.23*resolution_fraction[1]),
    //         threshold[2],
	// 		new RooRealVar("num1_254", "num", 25, 0.1, 100)
	// 	}
	// ));
	// 260
	info_map[1].insert(std::make_pair(
		"260",
		VoigtInfo{
			new RooRealVar("mean1_260", "mean", 26.2, 26.0, 26.4),
			new RooRealVar("gamma1_260", "gamma", 0.02, 0.0, 1.0),
			new RooConstVar("sigma1_260", "sigma", 0.25*resolution_fraction[1]),
            threshold[2],
			new RooRealVar("num1_260", "num", 25, 0.1, 200)
		}
	));
	// P.D.F.
	std::vector<RooAbsPdf*> pdf1;
	// extend P.D.F.
	std::vector<RooExtendPdf*> extend_pdf1;
	for (auto &[key, value] : info_map[1]) {
		// check
		if (!value.mean || !value.gamma || !value.sigma) {
			std::cerr << "Error: Missing variable for state 1, energy "
				<< key << "\n";
			return -1;
		}
		// construct P.D.F.
		pdf1.push_back(new AsymmetricVoigtian(
			("voigt1_"+key).c_str(), "voigt",
			x1,
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
	// background
	RooRealVar bkg1_root0("bkg1_root0", "root0", 17, 16, 18);
	RooRealVar bkg1_root1("bkg1_root1", "root1", 26, 25, 27);
	BackgroundPoly bkg1("bkg1", "bkg", x1, bkg1_root0, bkg1_root1);
	// background events' number
	RooRealVar num_bkg1("num_bkg1", "num", 10, 0, 30);
	RooExtendPdf ext_bkg1("ext_bkg1", "extend", bkg1, num_bkg1);
	model1_list.add(ext_bkg1);
	// model1
	RooAddPdf model1("model1", "model", model1_list);


	// 10Be 6 MeV
	RooDataHist hist2("dh2", "excited energy", RooArgList(x2), hist_ex+2);
	// RooDataSet set2("s2", "excited energy", fit_tree[2], RooArgSet(x));
	// construct P.D.F. information
	// 19.8
	info_map[2].insert(std::make_pair(
		"198",
		VoigtInfo {
			new RooRealVar("mean2_198", "mean", 19.8, 19.5, 20.0),
			new RooRealVar("gamma2_198", "gamma", 0.05, 0.0, 0.2),
			new RooConstVar("sigma2_198", "sigma", 0.11*resolution_fraction[2]),
			threshold[2],
			new RooRealVar("num2_198", "num", 15, 0.1, 100)
		}
	));
	// 20.1
	// info_map[2].insert(std::make_pair(
	// 	"201",
	// 	VoigtInfo {
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma2_201", "sigma", 0.11*resolution_fraction[2]),
	// 		threshold[2],
	// 		new RooRealVar("num2_201", "num", 15, 0.1, 100)
	// 	}
	// ));
	// 20.7
	info_map[2].insert(std::make_pair(
		"207",
		VoigtInfo {
			new RooRealVar("mean2_207", "mean", 20.7, 20.5, 20.8),
			new RooRealVar("gamma2_207", "gamma", 0.02, 0.0, 0.4),
			new RooConstVar("sigma2_207", "sigma", 0.12*resolution_fraction[2]),
			threshold[2],
			new RooRealVar("num2_207", "num", 50, 0.1, 150)
		}
	));
	// 20.9
	// info_map[2].insert(std::make_pair(
	// 	"209",
	// 	VoigtInfo {
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma2_209", "sigma", 0.13*resolution_fraction[2]),
	// 		threshold[2],
	// 		new RooRealVar("num2_209", "num", 40, 0.1, 150)
	// 	}
	// ));
	// 21.6
	info_map[2].insert(std::make_pair(
		"216",
		VoigtInfo {
			new RooRealVar("mean2_216", "mean", 21.6, 21.3, 21.8),
			new RooRealVar("gamma2_216", "gamma", 0.05, 0.0, 1.0),
			new RooConstVar("sigma2_216", "sigma", 0.14*resolution_fraction[2]),
			threshold[2],
			new RooRealVar("num2_216", "num", 50, 0.1, 150)
		}
	));
	// 22.3
	info_map[2].insert(std::make_pair(
		"223",
		VoigtInfo {
			new RooRealVar("mean2_223", "mean", 22.3, 22.0, 22.6),
			new RooRealVar("gamma2_223", "gamma", 0.04, 0.0, 1.0),
			new RooConstVar("sigma2_223", "sigma", 0.15*resolution_fraction[2]),
			threshold[2],
			new RooRealVar("num2_223", "num", 120, 0.1, 300)
		}
	));
	// 22.8
	// info_map[2].insert(std::make_pair(
	// 	"228",
	// 	VoigtInfo {
	// 		nullptr,
	// 		nullptr,
	// 		new RooConstVar("sigma2_228", "sigma", 0.17*resolution_fraction[2]),
	// 		threshold[2],
	// 		new RooRealVar("num2_228", "num", 10, 0.1, 200)
	// 	}
	// ));
	// 23.2
	info_map[2].insert(std::make_pair(
		"232",
		VoigtInfo {
			new RooRealVar("mean2_232", "mean", 23.2, 23.0, 23.5),
			new RooRealVar("gamma2_232", "gamma", 0.02, 0.0, 0.2),
			new RooConstVar("sigma2_232", "sigma", 0.17*resolution_fraction[2]),
			threshold[2],
			new RooRealVar("num2_232", "num", 50, 0.1, 200)
		}
	));
	// 23.6
	info_map[2].insert(std::make_pair(
		"236",
		VoigtInfo {
			new RooRealVar("mean2_236", "mean", 23.6, 23.3, 24.2),
			new RooRealVar("gamma2_236", "gamma", 0.05, 0.0, 0.2),
			new RooConstVar("sigma2_236", "sigma", 0.18*resolution_fraction[2]),
			threshold[2],
			new RooRealVar("num2_236", "num", 50, 0.1, 200)
		}
	));
	// 24.5
	info_map[2].insert(std::make_pair(
		"245",
		VoigtInfo {
			new RooRealVar("mean2_245", "mean", 24.5, 24.2, 24.9),
			new RooRealVar("gamma2_245", "gamma", 0.02, 0.0, 0.2),
			new RooConstVar("sigma2_245", "sigma", 0.19*resolution_fraction[2]),
			threshold[2],
			new RooRealVar("num2_245", "num", 50, 0.1, 200)
		}
	));
	// 25.4
	info_map[2].insert(std::make_pair(
		"254",
		VoigtInfo {
			new RooRealVar("mean2_254", "mean", 25.4, 25.2, 25.7),
			new RooRealVar("gamma2_254", "gamma", 0.02, 0.0, 1.0),
			new RooConstVar("sigma2_254", "sigma", 0.20*resolution_fraction[2]),
			threshold[2],
			new RooRealVar("num2_254", "num", 30, 0.1, 200)
		}
	));
	// 26.0
	info_map[2].insert(std::make_pair(
		"260",
		VoigtInfo {
			new RooRealVar("mean2_260", "mean", 26.2, 26.0, 26.4),
			new RooRealVar("gamma2_260", "gamma", 0.02, 0.0, 1.0),
			new RooConstVar("sigma2_260", "sigma", 0.21*resolution_fraction[2]),
			threshold[2],
			new RooRealVar("num2_260", "num", 30, 0.1, 200)
		}
	));

	// P.D.F.
	std::vector<RooAbsPdf*> pdf2;
	// extended P.D.F.
	std::vector<RooExtendPdf*> extend_pdf2;
	for (auto &[key, value] : info_map[2]) {
		// check
		if (!value.mean || !value.gamma || !value.sigma) {
			std::cerr << "Error: Missing variable for state 2, energy "
				<< key << "\n";
			return -1;
		}
		// construct P.D.F.
		pdf2.push_back(new AsymmetricVoigtian(
			("voigt2_"+key).c_str(), "voigt",
			x2,
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
	BackgroundPoly bkg2("bkg2", "bkg", x2, bkg2_root0, bkg2_root1);
	// background events' number
	RooRealVar num_bkg2("num_bkg2", "num", 10, 0, 30);
	RooExtendPdf ext_bkg2("ext_bkg2", "extend", bkg2, num_bkg2);
	model2_list.add(ext_bkg2);
	// model0
	RooAddPdf model2("model2", "model", model2_list);


	std::unique_ptr<RooFitResult> fit_result0{
		model0.fitTo(
			hist0,
			RooFit::PrintLevel(-1),
			RooFit::Save()
		)
	};
	fit_result0->Print();
	std::unique_ptr<RooFitResult> fit_result1{
		model1.fitTo(
			hist1,
			RooFit::PrintLevel(-1),
			RooFit::Save()
		)
	};
	fit_result1->Print();
	std::unique_ptr<RooFitResult> fit_result2{
		model2.fitTo(
			hist2,
			RooFit::PrintLevel(-1),
			RooFit::Save()
		)
	};
	fit_result2->Print();

	// show bar
	// std::vector<std::string> bar_names = {
	// 	"172", "185", "190", "193",
	// 	"201", "205", "216", "223",
	// 	"228", "232", "236", "245"
	// };
	// for (size_t i = 0; i < bar_names.size(); ++i) {
	// 	for (int state = 0; state < 3; ++state) {
	// 		auto search = info_map[state].find(bar_names[i]);
	// 		if (search == info_map[state].end()) {
	// 			std::cerr << "Warning: Could not find "
	// 				<< bar_names[i] << " in state " << state << "\n";
	// 			continue;
	// 		}
	// 		double mean = search->second.mean->getVal();
	// 		double event = search->second.num->getVal();
	// 		double relative_event = event / g_efficiency[state]->Eval(mean);
	// 		relative_bar[state].SetBinContent(i+1, relative_event);
	// 	}
	// }

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
		{0, "183"},
		{1, "183"},
		{0, "195"},
		{1, "195"},
		{0, "202"},
		{1, "207"},
		{2, "207"},
		{0, "216"},
		{1, "216"},
		{2, "216"},
		{0, "223"},
		{2, "223"},
		{0, "228"},
		{1, "228"},
		{0, "232"},
		{2, "232"},
		{0, "236"},
		{1, "236"},
		{2, "236"},
		{0, "245"},
		{1, "245"},
		{0, "254"},
		{2, "254"},
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
	RooPlot *frame0 = x0.frame();
	hist0.plotOn(
		frame0,
		RooFit::DrawOption("0"),
		RooFit::DataError(RooAbsData::None)
		// RooFit::FillColor(kCyan)
	);
	model0.plotOn(frame0, RooFit::LineColor(kRed));
	for (size_t i = 0; i < pdf0.size(); ++i) {
		PlotPdf(model0, frame0, *pdf0[i]);
	}
	model0.plotOn(
		frame0,
		RooFit::LineColor(kOrange),
		RooFit::Components(bkg0)
	);

	// frame 1
	RooPlot *frame1 = x1.frame();
	hist1.plotOn(
		frame1,
		RooFit::DrawOption("0"),
		RooFit::DataError(RooAbsData::None)
	);
	model1.plotOn(frame1, RooFit::LineColor(kRed));
	for (size_t i = 0; i < pdf1.size(); ++i) {
		PlotPdf(model1, frame1, *pdf1[i]);
	}
	model1.plotOn(
		frame1,
		RooFit::LineColor(kOrange),
		RooFit::Components(bkg1)
	);

	// frame 2
	RooPlot *frame2 = x2.frame();
	hist2.plotOn(
		frame2,
		RooFit::DrawOption("0"),
		RooFit::DataError(RooAbsData::None)
	);
	model2.plotOn(frame2, RooFit::LineColor(kRed));
	for (size_t i = 0; i < pdf2.size(); ++i) {
		PlotPdf(model2, frame2, *pdf2[i]);
	}
	model2.plotOn(
		frame2,
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
	opf.cd();
	for (int i = 0; i < 3; ++i) hist_ex[i].Write();
	frame0->Write("rh0");
	frame1->Write("rh1");
	frame2->Write("rh2");
	// for (size_t i = 0; i < 3; ++i) relative_bar[i].Write();
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
		TString::Format(
			"%s%sroofit-alone-result.png",
			kGenerateDataPath,
			kImageDir
		)
	);
	c1->Write("c1");


	// TCanvas *c2 = new TCanvas("c2", "", 1920, 1080);
	// gStyle->SetOptStat(0);
	// TAxis *xaxis = relative_bar[2].GetXaxis();
	// for (size_t i = 0; i < bar_names.size(); i++) {
	// 	// get same mean and gamma values
	// 	auto search = info_map[0].find(bar_names[i]);
	// 	if (search == info_map[0].end()) continue;
	// 	double mean = search->second.mean->getVal();
	// 	xaxis->SetNdivisions(13);
	// 	xaxis->ChangeLabel(
	// 		i+2, -1, -1, -1, -1, -1,
	// 		TString::Format("%.1lf", mean)
	// 	);
	// }
	// relative_bar[2].Draw("bar");
	// relative_bar[1].Draw("bar, same");
	// relative_bar[0].Draw("bar, same");
	// TLegend *legend = new TLegend(0.72, 0.75, 0.88, 0.85);
	// legend->AddEntry(relative_bar+0, "ground state", "f");
	// legend->AddEntry(relative_bar+1, "2+ state", "f");
	// legend->AddEntry(relative_bar+2, "~6MeV state", "f");
	// legend->Draw();
	// c2->SaveAs(
	// 	TString::Format("%s%srelative-events.png", kGenerateDataPath, kImageDir)
	// );
	// c2->Write("barc");

	opf.Close();

	return 0;
}