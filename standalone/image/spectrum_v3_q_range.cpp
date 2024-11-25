#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

#include "include/defs.h"

using namespace ribll;

constexpr int bins = 50;
constexpr double step = 0.01;


int FillSpectrumV2(
	std::vector<TH1F> &hist_ex0,
	std::vector<TH1F> &hist_ex1,
	std::vector<TH1F> &hist_ex2
) {
		// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody-10Be-dv2-2.root",
		kGenerateDataPath,
		kSpectrumDir
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
	int ppac_flag, taf_flag, bind, hole, target_flag;
	// PID in correct range? index 0 for origin, 1 for calculated
	int straight[2];
	// Q value
	double q[4];
	// excited energy ignore 10Be state
	double stateless_excited_energy[3][4];
	// setup input branches
	ipt->SetBranchAddress("ppac_flag", &ppac_flag);
	ipt->SetBranchAddress("taf_flag", &taf_flag);
	ipt->SetBranchAddress("bind", &bind);
	ipt->SetBranchAddress("hole", &hole);
	ipt->SetBranchAddress("target_flag", &target_flag);
	ipt->SetBranchAddress("straight", straight);
	ipt->SetBranchAddress("q", q);
	ipt->SetBranchAddress("stateless_excited_energy", stateless_excited_energy);


	// loop and fill events
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		// get data
		ipt->GetEntry(entry);
		// ignore bad events
		if (ppac_flag == 0) continue;
		if (taf_flag != 0) continue;
		if (target_flag != 1) continue;
		if (bind != 0) continue;
		if (hole != 0) continue;
		if (straight[0] != 3) continue;

		// ground state
		if (q[0] < -11 && q[0] > -13) {
			for (int i = 0; i < 30; ++i) {
				hist_ex0[i].Fill(stateless_excited_energy[0][0]);
			}
		}
		// first excited state
		if (q[0] < -14.5 && q[0] > -16) {
			for (int i = 0; i < 30; ++i) {
				hist_ex1[i].Fill(stateless_excited_energy[1][0]);
			}
		}
		// 6MeV state
		if (q[0] < -17 && q[0] > -20) {
			for (int i = 0; i < 30; ++i) {
				hist_ex2[i].Fill(stateless_excited_energy[2][0]);
			}
		}
	}
	// close file
	ipf.Close();
	return 0;
}



int FillSpectrumV3(
	std::vector<TH1F> &hist_ex0,
	std::vector<TH1F> &hist_ex1,
	std::vector<TH1F> &hist_ex2
) {
		// input file name
	TString input_file_name = TString::Format(
		"%s%sC14-10Be-4He-2H-v3.root",
		kGenerateDataPath,
		kSpectrumDir
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
	// Q value and excited energy
	double q, excited_energy;
	// setup input branches
	ipt->SetBranchAddress("valid", &valid);
	ipt->SetBranchAddress("taf_flag", &taf_flag);
	ipt->SetBranchAddress("be_state", &be_state);
	ipt->SetBranchAddress("q", &q);
	ipt->SetBranchAddress("excited_energy_target", &excited_energy);

	// loop and fill events
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		// get data
		ipt->GetEntry(entry);
		// ignore bad events
		if (valid != 0) continue;
		if (taf_flag == 2) continue;

		for (size_t i = 0; i < 30; ++i) {
			if (be_state == 0) {
				hist_ex0[i].Fill(excited_energy);
			} else if (be_state == 1) {
				hist_ex1[i].Fill(excited_energy);
			} else if (be_state == 2) {
				hist_ex2[i].Fill(excited_energy);
			}
		}
	}
	// close file
	ipf.Close();
	return 0;
}


int main() {
	// output file name
	TString output_file_name = TString::Format(
		"%s%sspectrum-v3.root",
		kGenerateDataPath,
		kImageDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// ground state excited energy spectrum, change Q range and start point
	std::vector<TH1F> hist_ex0_qs;
	std::vector<TH1F> hist_ex1_qs;
	std::vector<TH1F> hist_ex2_qs;

	// set histograms parameters, changes start point
	// ground state
	for (size_t i = 0; i < 30; ++i) {
		hist_ex0_qs.emplace_back(
			TString::Format("hex0r%ld", i),
			TString::Format(
				"ex0 %.2lf to %.2lf",
				12.0 + step*i,
				27.0 + step*i
			),
			bins, 12.0 + step*i, 27.0 + step*i
		);
	}
	// first excited state
	for (size_t i = 0; i < 30; ++i) {
		hist_ex1_qs.emplace_back(
			TString::Format("hex1r%ld", i),
			TString::Format(
				"ex1 %.2lf to %.2lf",
				12.0 + step*i,
				27.0 + step*i
			),
			bins, 12.0+step*i, 27.0+step*i
		);
	}
	// 6MeV state
	for (size_t i = 0; i < 30; ++i) {
		hist_ex2_qs.emplace_back(
			TString::Format("hex2r%ld", i),
			TString::Format(
				"ex2 %.2lf to %.2lf",
				12.0 + step*i,
				27.0 + step*i
			),
			bins, 12.0+step*i, 27.0+step*i
		);
	}

	if (FillSpectrumV2(
		hist_ex0_qs, hist_ex1_qs, hist_ex2_qs
	)) {
		std::cerr << "Error: Fill spectrum V2 failed.\n";
		return -1;
	}


	if (FillSpectrumV3(
		hist_ex0_qs, hist_ex1_qs, hist_ex2_qs
	)) {
		std::cerr << "Error: Fill spectrum V3 failed.\n";
		return -1;
	}

	// save images
	opf.cd();
	for (size_t i = 0; i < 30; ++i) {
		hist_ex0_qs[i].Write();
	}
	for (size_t i = 0; i < 30; ++i) {
		hist_ex1_qs[i].Write();
	}
	for (size_t i = 0; i < 30; ++i) {
		hist_ex2_qs[i].Write();
	}

	// save as pdf file, ground state
	TCanvas *c1 = new TCanvas("c1", "c1", 1920, 1080);
	c1->cd();
	TString pdf_file_name_0 = TString::Format(
		"%s%sspectrum3-ex0-qs.pdf",
		kGenerateDataPath,
		kImageDir
	);
	c1->Print(pdf_file_name_0+"[");
	for (size_t i = 0; i < 30; ++i) {
		hist_ex0_qs[i].Draw();
		c1->Print(pdf_file_name_0);
	}
	c1->Print(pdf_file_name_0+"]");
	// first excited state
	c1->cd();
	TString pdf_file_name_1 = TString::Format(
		"%s%sspectrum3-ex1-qs.pdf",
		kGenerateDataPath,
		kImageDir
	);
	c1->Print(pdf_file_name_1+"[");
	for (size_t i = 0; i < 30; ++i) {
		hist_ex1_qs[i].Draw();
		c1->Print(pdf_file_name_1);
	}
	c1->Print(pdf_file_name_1+"]");
	// 6MeV state
	c1->cd();
	TString pdf_file_name_2 = TString::Format(
		"%s%sspectrum3-ex2-qs.pdf",
		kGenerateDataPath,
		kImageDir
	);
	c1->Print(pdf_file_name_2+"[");
	for (size_t i = 0; i < 30; ++i) {
		hist_ex2_qs[i].Draw();
		c1->Print(pdf_file_name_2);
	}
	c1->Print(pdf_file_name_2+"]");

	// close files
	opf.Close();
	return 0;
}