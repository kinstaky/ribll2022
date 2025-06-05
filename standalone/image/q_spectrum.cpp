/*
 * 绘制拟合最终的 Q 值谱
 */

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

int FillV2(TH1F &hq, TH1F &hq0, TH1F &hq1, TH1F &hq2) {
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
	double c_kinetic[4], q[4];
	// setup input branches
	ipt->SetBranchAddress("ppac_flag", &ppac_flag);
	ipt->SetBranchAddress("taf_flag", &taf_flag);
	ipt->SetBranchAddress("target_flag", &target_flag);
	ipt->SetBranchAddress("bind", &bind);
	ipt->SetBranchAddress("hole", &hole);
	ipt->SetBranchAddress("straight", straight);
	ipt->SetBranchAddress("c_kinetic", c_kinetic);
	ipt->SetBranchAddress("q", q);

	// loop
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);
		if (ppac_flag  == 0) continue;
		if (taf_flag != 0) continue;
		if ((target_flag & 1) == 0) continue;
		if (bind != 0) continue;
		if (hole > 0) continue;
		// straight PID
		if (straight[0] != 3) continue;
		// beam kinetic energy
		if (c_kinetic[0] < 360.0) continue;

		hq.Fill(q[0]);
		hq0.Fill(q[0]);
		hq1.Fill(q[0]);
		hq2.Fill(q[0]);
	}

	// close file
	ipf.Close();
	return 0;
}


int FillV3(TH1F &hq, TH1F &hq0, TH1F &hq1, TH1F &hq2) {
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
	int valid, taf_flag;
	// Q value
	double q;
	// setup input branches
	ipt->SetBranchAddress("valid", &valid);
	ipt->SetBranchAddress("taf_flag", &taf_flag);
	ipt->SetBranchAddress("q", &q);

	// loop
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		if (valid != 0) continue;
		// if (taf_flag != 1) continue;
		hq.Fill(q);
		hq0.Fill(q);
		hq1.Fill(q);
		hq2.Fill(q);
	}

	ipf.Close();
	return 0;
}


inline double Gaus(double x, double a, double mean, double sigma) {
	return a * exp(-0.5 * pow((x - mean) / sigma, 2.0));
}

double FitQ(double *x, double *par) {
	return Gaus(x[0], par[0], par[3], par[4])
		+ Gaus(x[0], par[1], par[3]-3.368, par[5])
		+ Gaus(x[0], par[2], par[3]-6.179, par[6]);
}


int main() {
	// output file name
	TString output_file_name = TString::Format(
		"%s%sq-final.root",
		kGenerateDataPath,
		kSpectrumDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// Q value spectrum
	TH1F hq("hq", "Q value", 90, -23, -8);
	TH1F hq0("hq0", "Q value for state 0", 12, -13.0, -11.0);
	TH1F hq1("hq1", "Q value for state 1", 6, -16.0, -15.0);
	TH1F hq2("hq2", "Q value for state 2", 12, -19.5, -17.5);
	// output data
	QEvent event;
	// output tree
	TTree opt("tree", "Q");
	// setup output branches
	opt.Branch("version", &event.version, "ver/I");
	opt.Branch("q", &event.q, "q/D");

	// if (FillV2(hq, hq0, hq1, hq2)) return -1;
	if (FillV3(hq, hq0, hq1, hq2)) return -1;

	// result is -12.5247, -15.5002, -17.9949
	constexpr double q_init_param[9] = {
		80, -12, 1,
		200, -15.5, 0.6,
		100, -18, 0.9
	};
	// fit Q spectrum
	TF1 *fq = new TF1("fq", "gaus(0)+gaus(3)+gaus(6)", -23, -8);
	fq->SetParameters(q_init_param);
	fq->SetParLimits(1, -14, -11);
	fq->SetParLimits(4, -17, -14);
	fq->SetParLimits(7, -19, -17);
	fq->SetNpx(1000);
	hq.Fit(fq, "R+");
	std::cout << fq->GetParameter(1) << ", "
		<< fq->GetParameter(4) << ", "
		<< fq->GetParameter(7) << "\n";
	for (size_t i = 0; i < 3; ++i) {
		TF1 *f1 = new TF1(TString::Format("fq%ld", i), "gaus", -23, -8);
		f1->SetLineColor(kBlue);
		f1->SetLineWidth(3);
		f1->SetLineStyle(2);
		f1->SetNpx(1000);
		f1->SetParameters(fq->GetParameters()+3*i);
		hq.GetListOfFunctions()->Add(f1);
	}

	// // result is -12.1234
	// constexpr double q_init_param[7] = {
	// 	80, 200, 100,
	// 	-12,
	// 	1, 0.6, 0.9,
	// };
	// // fit Q spectrum
	// TF1 *fq = new TF1("fq", FitQ, -23, -8, 7);
	// fq->SetParameters(q_init_param);
	// fq->SetParLimits(3, -13, -11);
	// fq->SetNpx(1000);
	// hq.Fit(fq, "R+");
	// std::cout << fq->GetParameter(3) << "\n";
	// for (size_t i = 0; i < 3; ++i) {
	// 	TF1 *f1 = new TF1(TString::Format("fq%ld", i), "gaus", -23, -8);
	// 	f1->SetLineColor(kBlue);
	//	f1->SetNpx(1000);
	// 	f1->SetParameter(0, fq->GetParameter(0+i));
	// 	if (i == 0) {
	// 		f1->SetParameter(1, fq->GetParameter(3));
	// 	} else if (i == 1) {
	// 		f1->SetParameter(1, fq->GetParameter(3)-3.368);
	// 	} else if (i == 2) {
	// 		f1->SetParameter(1, fq->GetParameter(3)-6.179);
	// 	}
	// 	f1->SetParameter(2, fq->GetParameter(4+i));
	// 	f1->SetLineWidth(3);
	// 	f1->SetLineStyle(2);
	// 	hq.GetListOfFunctions()->Add(f1);
	// }

	TCanvas *c1 = new TCanvas("c1", "c1", 1520, 855);
	c1->cd();
	gStyle->SetOptTitle(false);
	hq.SetStats(false);
	hq.GetXaxis()->SetLabelFont(132);
	// hq.GetXaxis()->SetMaxDigits(3);
	hq.GetXaxis()->SetLabelSize(0.06);
	hq.GetYaxis()->SetLabelFont(132);
	// hq.GetYaxis()->SetMaxDigits(3);
	hq.GetYaxis()->SetLabelSize(0.06);
	hq.Draw();
	hq0.SetFillColor(kBlack);
	hq0.SetFillStyle(3354);
	hq0.Draw("same");
	hq1.SetFillColor(kBlack);
	hq1.SetFillStyle(3354);
	hq1.Draw("same");
	hq2.SetFillColor(kBlack);
	hq2.SetFillStyle(3354);
	hq2.Draw("same");
	c1->Print(TString::Format(
		"%s%sq-final.png", kGenerateDataPath, kImageDir
	));

	// save
	opf.cd();
	hq.Write();
	hq0.Write();
	hq1.Write();
	hq2.Write();
	c1->Write();
	// close files
	opf.Close();
	return 0;
}