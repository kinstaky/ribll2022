#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "include/defs.h"

using namespace ribll;

struct QEvent {
	int version;
	double q;
};

int FillV2(
	std::vector<TH1F> &hex0,
	std::vector<TH1F> &hex1,
	std::vector<TH1F> &hex2
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
	int valid;
	double q[4], stateless_excited_energy[3][4];
	// setup input branches
	ipt->SetBranchAddress("valid", &valid);
	ipt->SetBranchAddress("q", q);
	ipt->SetBranchAddress(
		"stateless_excited_energy", stateless_excited_energy
	);

	// loop
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);
		if (valid) continue;

		int index = -1;
		double ex = -1.0;
		if (q[0] < -11 && q[0] > -13) {
			ex = stateless_excited_energy[0][0];
			index = 0;
		} else if (q[0] < -14.5 && q[0] > -16) {
			ex = stateless_excited_energy[1][0];
			index = 1;
		} else if (q[0] < -17 && q[0] > -20) {
			ex = stateless_excited_energy[2][0];
			index = 2;
		}
		if (ex > 0) {
			if (index == 0) {
				for (size_t i = 0; i < hex0.size(); ++i) {
					hex0[i].Fill(ex);
				}
			} else if (index == 1) {
				for (size_t i = 0; i < hex1.size(); ++i) {
					hex1[i].Fill(ex);
				}
			}
			if (index == 2) {
				for (size_t i = 0; i < hex2.size(); ++i) {
					hex2[i].Fill(ex);
				}
			}
		}
	}

	// close file
	ipf.Close();
	return 0;
}


int FillV3(
	std::vector<TH1F> &hex0,
	std::vector<TH1F> &hex1,
	std::vector<TH1F> &hex2
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
	int valid, taf_flag, be_state;
	// Q value
	double excited_energy;
	// setup input branches
	ipt->SetBranchAddress("valid", &valid);
	ipt->SetBranchAddress("taf_flag", &taf_flag);
	ipt->SetBranchAddress("be_state", &be_state);
	ipt->SetBranchAddress("excited_energy_target", &excited_energy);

	// loop
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		if (valid != 0) continue;
		// if (taf_flag != 1) continue;
		if (be_state == 0) {
			for (size_t i = 0; i < hex0.size(); ++i) {
				hex0[i].Fill(excited_energy);
			}
		} else if (be_state == 1) {
			for (size_t i = 0; i < hex1.size(); ++i) {
				hex1[i].Fill(excited_energy);
			}
		} else if (be_state == 2) {
			for (size_t i = 0; i < hex2.size(); ++i) {
				hex2[i].Fill(excited_energy);
			}
		}
	}

	ipf.Close();
	return 0;
}


int main() {
	// output file name
	TString output_file_name = TString::Format(
		"%s%sfinal-excitation-spectrum-bin.root",
		kGenerateDataPath,
		kSpectrumDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// excitation spectrum
	std::vector<TH1F> hist_ex[3];
	for (size_t i = 0; i < 3; ++i) {
		for (size_t j = 0; j < 30; ++j) {
			hist_ex[i].emplace_back(
                Form("hex%ld_%ld", i, j),
                Form(
					"excitation spectrum %ld, from %.2f to %.2f",
					i, 12.0+0.01*j, 27.0+0.01*j
				),
                60, 12.0+0.01*j, 27.0+0.01*j
            );
		}
	}

	// if (FillV2(hist_ex[0], hist_ex[1], hist_ex[2])) return -1;
	if (FillV3(hist_ex[0], hist_ex[1], hist_ex[2])) return -1;

	TCanvas *c1 = new TCanvas("c1", "c1", 1920, 1080);
	c1->cd();
	TString pdf_file_name = TString::Format(
		"%s%sspectrum-final-bin.pdf",
		kGenerateDataPath,
		kImageDir
	);
	c1->Print(pdf_file_name+"[");
	for (size_t i = 0; i < 3; ++i) {
		for (auto &hist : hist_ex[i]) {
			hist.Draw();
			c1->Print(pdf_file_name);
		}
	}
	c1->Print(pdf_file_name+"]");

	// save
	opf.cd();
	for (size_t i = 0; i < 3; ++i) {
		for (auto &hist : hist_ex[i]) hist.Write();
	}
	// close files
	opf.Close();
	return 0;
}