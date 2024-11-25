#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLine.h>

#include "include/defs.h"

using namespace ribll;

constexpr int bins = 50;
constexpr double step = 0.01;
constexpr size_t ex0_count = 7;
constexpr size_t ex1_count = 8;
constexpr size_t ex2_count = 12;

int main(int argc, char **argv) {
	std::string suffix = "";
	if (argc > 1) {
		if (argv[1][0] == '-' && argv[1][1] == 'h') {
			std::cout << "Usage: " << argv[0] << " [suffix]\n"
				<< "  suffix        ROOT file suffix tag\n";
			return -1;
		}
		suffix = std::string(argv[1]);
	}

	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody%s-2.root",
		kGenerateDataPath,
		kSpectrumDir,
		suffix.empty() ? "" : ("-"+suffix).c_str()
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

	// output file name
	TString output_file_name = TString::Format(
		"%s%sspectrum2%s.root",
		kGenerateDataPath,
		kShowDir,
		suffix.empty() ? "" : ("-"+suffix).c_str()
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// ground state excited energy spectrum, change Q range and start point
	TH1F hist_ex0_qs[ex0_count][32];
	TH1F hist_ex1_qs[ex1_count][32];
	TH1F hist_ex2_qs[ex2_count][32];
	std::vector<TH1F> hist_ex0_q, hist_ex1_q, hist_ex2_q;

	// set histograms parameters, changes Q range and start point
	// ground state
	std::pair<double, double> ex0_qs_range[ex0_count] = {
		{-13.0, -11.0},
		{-13.0, -11.2},
		{-13.0, -10.8},
		{-13.0, -10.0},
		{-13.5, -11.0},
		{-14.0, -11.0},
		{-14.0, -10.0}
	};
	for (size_t i = 0; i < ex0_count; ++i) {
		for (int j = 0; j < 30; ++j) {
			hist_ex0_qs[i][j] = TH1F(
				TString::Format("hex0q%lde%d", i, j),
				TString::Format(
					"ground state Q %.2lf to %.2lf, ex %.2lf to %.2lf",
					ex0_qs_range[i].first,
					ex0_qs_range[i].second,
					12.0 + step*j,
					27.0 + step*j
				),
				bins, 12.0 + step*j, 27.0 + step*j
			);
		}
	}
	// first excited state
	std::pair<double, double> ex1_qs_range[ex1_count] = {
		{-16.5, -14.5},
		{-16.0, -14.5},
		{-15.5, -14.5},
		{-16.2, -14.5},
		{-15.7, -14.5},
		{-16.5, -15.0},
		{-16.0, -15.0},
		{-15.5, -15.0}
	};
	for (size_t i = 0; i < ex1_count; ++i) {
		for (int j = 0; j < 30; ++j) {
			hist_ex1_qs[i][j] = TH1F(
				TString::Format("hex1q%lde%d", i, j),
				TString::Format(
					"2+ state Q %.2lf to %.2lf, ex %.2lf to %.2lf",
					ex1_qs_range[i].first,
					ex1_qs_range[i].second,
					12.0 + step*j,
					27.0 + step*j
				),
				bins, 12.0+step*j, 27.0+step*j
			);
		}
	}
	// 6MeV state
	std::pair<double, double> ex2_qs_range[ex2_count] = {
		{-20.0, -17.0},
		{-19.5, -17.0},
		{-19.0, -17.0},
		{-18.5, -17.0},
		{-20.0, -17.5},
		{-19.5, -17.5},
		{-19.0, -17.5},
		{-18.5, -17.5},
		{-20.0, -18.0},
		{-19.5, -18.0},
		{-19.0, -18.0},
		{-18.5, -18.0}
	};
	for (size_t i = 0; i < ex2_count; ++i) {
		for (int j = 0; j < 30; ++j) {
			hist_ex2_qs[i][j] = TH1F(
				TString::Format("hex2q%lde%d", i, j),
				TString::Format(
					"6MeV state Q %.2lf to %.2lf, ex %.2lf to %.2lf",
					ex2_qs_range[i].first,
					ex2_qs_range[i].second,
					12.0 + step*j,
					27.0 + step*j
				),
				bins, 12.0+step*j, 27.0+step*j
			);
		}
	}

	// excited energy from different Q
	std::vector<std::pair<double, double>> ex0_q_range;
	ex0_q_range.push_back(std::make_pair(-13.0, -11.0));	// reference
	ex0_q_range.push_back(std::make_pair(-13.5, -11.0));
	ex0_q_range.push_back(std::make_pair(-13.0, -11.0));
	ex0_q_range.push_back(std::make_pair(-12.5, -11.0));
	ex0_q_range.push_back(std::make_pair(-12.0, -11.0));
	ex0_q_range.push_back(std::make_pair(-11.5, -11.0));
	ex0_q_range.push_back(std::make_pair(-13.5, -11.5));
	ex0_q_range.push_back(std::make_pair(-13.0, -11.5));
	ex0_q_range.push_back(std::make_pair(-12.5, -11.5));
	ex0_q_range.push_back(std::make_pair(-12.0, -11.5));
	ex0_q_range.push_back(std::make_pair(-13.5, -10.5));
	ex0_q_range.push_back(std::make_pair(-13.0, -10.5));
	ex0_q_range.push_back(std::make_pair(-12.5, -10.5));
	ex0_q_range.push_back(std::make_pair(-12.0, -10.5));
	ex0_q_range.push_back(std::make_pair(-11.5, -10.5));
	for (size_t i = 0; i < ex0_q_range.size(); ++i) {
		hist_ex0_q.emplace_back(
			TString::Format("hex0q%ld", i),
			TString::Format(
				"ex0 Q %.2lf to %.2lf",
				ex0_q_range[i].first,
				ex0_q_range[i].second
			),
			50, 12.1, 27.1
		);
	}


	std::vector<std::pair<double, double>> ex1_q_range;
	ex1_q_range.push_back(std::make_pair(-16.0, -14.5));	// reference
	ex1_q_range.push_back(std::make_pair(-17.0, -14.5));
	ex1_q_range.push_back(std::make_pair(-16.5, -14.5));
	ex1_q_range.push_back(std::make_pair(-16.0, -14.5));
	ex1_q_range.push_back(std::make_pair(-15.7, -14.5));
	ex1_q_range.push_back(std::make_pair(-15.5, -14.5));
	ex1_q_range.push_back(std::make_pair(-15.2, -14.5));
	ex1_q_range.push_back(std::make_pair(-15.0, -14.5));
	ex1_q_range.push_back(std::make_pair(-17.0, -15.0));
	ex1_q_range.push_back(std::make_pair(-16.5, -15.0));
	ex1_q_range.push_back(std::make_pair(-16.0, -15.0));
	ex1_q_range.push_back(std::make_pair(-15.7, -15.0));
	ex1_q_range.push_back(std::make_pair(-15.5, -15.0));
	ex1_q_range.push_back(std::make_pair(-17.0, -14.0));
	ex1_q_range.push_back(std::make_pair(-16.5, -14.0));
	ex1_q_range.push_back(std::make_pair(-16.0, -14.0));
	ex1_q_range.push_back(std::make_pair(-15.7, -14.0));
	ex1_q_range.push_back(std::make_pair(-15.5, -14.0));
	ex1_q_range.push_back(std::make_pair(-15.2, -14.0));
	ex1_q_range.push_back(std::make_pair(-15.0, -14.0));
	ex1_q_range.push_back(std::make_pair(-14.5, -14.0));
	for (size_t i = 0; i < ex1_q_range.size(); ++i) {
		hist_ex1_q.emplace_back(
			TString::Format("hex1q%ld", i),
			TString::Format(
				"ex1 Q %.2lf to %.2lf",
				ex1_q_range[i].first,
				ex1_q_range[i].second
			),
			50, 12.1, 27.1
		);
	}

	std::vector<std::pair<double, double>> ex2_q_range;
	ex2_q_range.push_back(std::make_pair(-20.0, -17.0));    // reference
	ex2_q_range.push_back(std::make_pair(-20.5, -17.0));
	ex2_q_range.push_back(std::make_pair(-20.0, -17.0));
	ex2_q_range.push_back(std::make_pair(-19.5, -17.0));
	ex2_q_range.push_back(std::make_pair(-19.0, -17.0));
	ex2_q_range.push_back(std::make_pair(-18.5, -17.0));
	ex2_q_range.push_back(std::make_pair(-18.0, -17.0));
	ex2_q_range.push_back(std::make_pair(-17.5, -17.0));
	ex2_q_range.push_back(std::make_pair(-20.5, -16.5));
	ex2_q_range.push_back(std::make_pair(-20.0, -16.5));
	ex2_q_range.push_back(std::make_pair(-19.5, -16.5));
	ex2_q_range.push_back(std::make_pair(-19.0, -16.5));
	ex2_q_range.push_back(std::make_pair(-18.5, -16.5));
	ex2_q_range.push_back(std::make_pair(-18.0, -16.5));
	ex2_q_range.push_back(std::make_pair(-17.5, -16.5));
	ex2_q_range.push_back(std::make_pair(-17.0, -16.5));
	ex2_q_range.push_back(std::make_pair(-20.5, -17.5));
	ex2_q_range.push_back(std::make_pair(-20.0, -17.5));
	ex2_q_range.push_back(std::make_pair(-19.5, -17.5));
	ex2_q_range.push_back(std::make_pair(-19.0, -17.5));
	ex2_q_range.push_back(std::make_pair(-18.5, -17.5));
	ex2_q_range.push_back(std::make_pair(-18.0, -17.5));
	ex2_q_range.push_back(std::make_pair(-20.0, -18.0));
	ex2_q_range.push_back(std::make_pair(-19.5, -18.0));
	ex2_q_range.push_back(std::make_pair(-19.0, -18.0));
	ex2_q_range.push_back(std::make_pair(-18.5, -18.0));
	ex2_q_range.push_back(std::make_pair(-20.0, -18.5));
	ex2_q_range.push_back(std::make_pair(-19.5, -18.5));
	ex2_q_range.push_back(std::make_pair(-19.0, -18.5));
	for (size_t i = 0; i < ex2_q_range.size(); ++i) {
		hist_ex2_q.emplace_back(
			TString::Format("hex2q%ld", i),
			TString::Format(
				"ex2 Q %.2lf to %.2lf",
				ex2_q_range[i].first,
				ex2_q_range[i].second
			),
			50, 12.1, 27.1
		);
	}


	// loop and fill events
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		// get data
		ipt->GetEntry(entry);
		// ignore bad events
		if ((ppac_flag & 1) == 0) continue;
		if (taf_flag != 0) continue;
		if (target_flag != 1) continue;
		if (bind != 0) continue;
		if (hole != 0) continue;
		if (straight[0] != 3) continue;

		// ground state
		for (size_t i = 0; i < ex0_count; ++i) {
			if (q[0] < ex0_qs_range[i].first) continue;
			if (q[0] > ex0_qs_range[i].second) continue;
			for (int j = 0; j < 30; ++j) {
				hist_ex0_qs[i][j].Fill(stateless_excited_energy[0][0]);
			}
		}
		// first excited state
		for (size_t i = 0; i < ex1_count; ++i) {
			if (q[0] < ex1_qs_range[i].first) continue;
			if (q[0] > ex1_qs_range[i].second) continue;
			for (int j = 0; j < 30; ++j) {
				hist_ex1_qs[i][j].Fill(stateless_excited_energy[1][0]);
			}
		}
		// 6MeV state
		for (size_t i = 0; i < ex2_count; ++i) {
			if (q[0] < ex2_qs_range[i].first) continue;
			if (q[0] > ex2_qs_range[i].second) continue;
			for (int j = 0; j < 30; ++j) {
				hist_ex2_qs[i][j].Fill(stateless_excited_energy[2][0]);
			}
		}

		// ground state
		for (size_t i = 0; i < ex0_q_range.size(); ++i) {
			if (q[0] > ex0_q_range[i].first && q[0] < ex0_q_range[i].second) {
				hist_ex0_q[i].Fill(stateless_excited_energy[0][0]);
			}
		}
		// first excited state
		for (size_t i = 0; i < ex1_q_range.size(); ++i) {
			if (q[0] > ex1_q_range[i].first && q[0] < ex1_q_range[i].second) {
				hist_ex1_q[i].Fill(stateless_excited_energy[1][0]);
			}
		}
		// 6MeV state
		for (size_t i = 0; i < ex2_q_range.size(); ++i) {
			if (q[0] > ex2_q_range[i].first && q[0] < ex2_q_range[i].second) {
				hist_ex2_q[i].Fill(stateless_excited_energy[2][0]);
			}
		}
	}

	// save images
	for (size_t i = 0; i < ex0_count; ++i) {
		for (int j = 0; j < 30; ++j) {
			hist_ex0_qs[i][j].Write();
		}
	}
	for (size_t i = 0; i < ex1_count; ++i) {
		for (int j = 0; j < 30; ++j) {
			hist_ex1_qs[i][j].Write();
		}
	}
	for (size_t i = 0; i < ex2_count; ++i) {
		for (int j = 0; j < 30; ++j) {
			hist_ex2_qs[i][j].Write();
		}
	}

	// save as pdf file, ground state
	TCanvas *c1 = new TCanvas("c1", "c1", 1920, 1080);
	c1->cd();

	// TString pdf_file_name_0 = TString::Format(
	// 	"%simage/show-spectrum2%s-ex0-qs.pdf",
	// 	kGenerateDataPath,
	// 	suffix.empty() ? "" : ("-"+suffix).c_str()
	// );
	// c1->Print(pdf_file_name_0+"[");
	// for (size_t i = 0; i < ex0_count; ++i) {
	// 	for (int j = 0; j < 30; ++j) {
	// 		hist_ex0_qs[i][j].Draw();
	// 		c1->Print(pdf_file_name_0);
	// 	}
	// }
	// c1->Print(pdf_file_name_0+"]");
	// // first excited state
	// c1->cd();
	// TString pdf_file_name_1 = TString::Format(
	// 	"%simage/show-spectrum2%s-ex1-qs.pdf",
	// 	kGenerateDataPath,
	// 	suffix.empty() ? "" : ("-"+suffix).c_str()
	// );
	// c1->Print(pdf_file_name_1+"[");
	// for (size_t i = 0; i < ex1_count; ++i) {
	// 	for (int j = 0; j < 30; ++j) {
	// 		hist_ex1_qs[i][j].Draw();
	// 		c1->Print(pdf_file_name_1);
	// 	}
	// }
	// c1->Print(pdf_file_name_1+"]");
	// // 6MeV state
	// c1->cd();
	// TString pdf_file_name_2 = TString::Format(
	// 	"%simage/show-spectrum2%s-ex2-qs.pdf",
	// 	kGenerateDataPath,
	// 	suffix.empty() ? "" : ("-"+suffix).c_str()
	// );
	// c1->Print(pdf_file_name_2+"[");
	// for (size_t i = 0; i < ex2_count; ++i) {
	// 	for (int j = 0; j < 30; ++j) {
	// 		hist_ex2_qs[i][j].Draw();
	// 		c1->Print(pdf_file_name_2);
	// 	}
	// }
	// c1->Print(pdf_file_name_2+"]");

	// use THStack for better visualization
	std::vector<THStack> ex0_q_stack;
	for (size_t i = 1; i < ex0_q_range.size(); ++i) {
		ex0_q_stack.emplace_back(
			TString::Format("ex0s%ld", i),
			hist_ex0_q[i].GetTitle()
		);
		ex0_q_stack[i-1].Add(&(hist_ex0_q[i]));
		ex0_q_stack[i-1].Add(&(hist_ex0_q[0]));
	}
	// lines representing peaks
	std::vector<double> ex0_peaks = {
		13.7, 14.82, 15.59, 16.5, 17.39, 18.23, 19.28, 19.78,
		20.56, 21.23, 21.67, 22.72, 23.91
	};
	std::vector<TLine> ex0_lines;
	for (const auto &peak : ex0_peaks) {
		ex0_lines.emplace_back(
			peak, 0.0, peak, 20.0
		);
		ex0_lines.back().SetLineStyle(2);
	}
	// ground state
	TString pdf_file_name_0_q = TString::Format(
		"%simage/show-spectrum2%s-ex0-q.pdf",
		kGenerateDataPath,
		suffix.empty() ? "" : ("-"+suffix).c_str()
	);
	c1->Print(pdf_file_name_0_q+"[");
	hist_ex0_q[0].SetLineColor(kRed);
	for (auto &stack : ex0_q_stack) {
		stack.Draw("nostack");
		for (auto &line : ex0_lines) line.Draw();
		c1->Print(pdf_file_name_0_q);
	}
	for (size_t i = 1; i < hist_ex0_q.size(); ++i) {
		hist_ex0_q[i].Draw();
		c1->Print(pdf_file_name_0_q);
	}
	c1->Print(pdf_file_name_0_q+"]");

	std::vector<THStack> ex1_q_stack;
	for (size_t i = 1; i < hist_ex1_q.size(); ++i) {
		ex1_q_stack.emplace_back(
			TString::Format("ex1s%ld", i),
			hist_ex1_q[i].GetTitle()
		);
		ex1_q_stack.back().Add(&(hist_ex1_q[i]));
		ex1_q_stack.back().Add(&(hist_ex1_q[0]));
	}
	// lines representing peaks
	std::vector<double> ex1_peaks = {
		17.39, 18.23, 18.71, 19.28, 19.78,
		20.56, 21.23, 21.67, 22.28, 22.72,
		23.18, 23.91, 24.65
	};
	std::vector<TLine> ex1_lines;
	for (const auto &peak : ex1_peaks) {
		ex1_lines.emplace_back(
			peak, 0.0, peak, 40.0
		);
		ex1_lines.back().SetLineStyle(2);
	}
	// first excited state
	TString pdf_file_name_1_q = TString::Format(
		"%simage/show-spectrum2%s-ex1-q.pdf",
		kGenerateDataPath,
		suffix.empty() ? "" : ("-"+suffix).c_str()
	);
	c1->Print(pdf_file_name_1_q+"[");
	hist_ex1_q[0].SetLineColor(kRed);
	for (auto &stack : ex1_q_stack) {
		stack.Draw("nostack");
		for (auto &line : ex1_lines) line.Draw();
		c1->Print(pdf_file_name_1_q);
	}
	for (size_t i = 1; i < hist_ex1_q.size(); ++i) {
		hist_ex1_q[i].Draw();
		c1->Print(pdf_file_name_1_q);
	}
	c1->Print(pdf_file_name_1_q+"]");

	std::vector<THStack> ex2_q_stack;
	for (size_t i = 1; i < hist_ex2_q.size(); ++i) {
		ex2_q_stack.emplace_back(
			TString::Format("ex2s%ld", i),
			hist_ex2_q[i].GetTitle()
		);
		ex2_q_stack.back().Add(&(hist_ex2_q[i]));
		ex2_q_stack.back().Add(&(hist_ex2_q[0]));
	}
	// lines representing peaks
	std::vector<double> ex2_peaks = {
		19.78, 20.56, 21.23, 21.67, 22.28,
		23.18, 23.91, 24.65
	};
	std::vector<TLine> ex2_lines;
	for (const auto &peak : ex2_peaks) {
		ex2_lines.emplace_back(
			peak, 0.0, peak, 50.0
		);
		ex2_lines.back().SetLineStyle(2);
	}
	// 6MeV state
	TString pdf_file_name_2_q = TString::Format(
		"%simage/show-spectrum2%s-ex2-q.pdf",
		kGenerateDataPath,
		suffix.empty() ? "" : ("-"+suffix).c_str()
	);
	c1->Print(pdf_file_name_2_q+"[");
	hist_ex2_q[0].SetLineColor(kRed);
	for (auto &stack : ex2_q_stack) {
		stack.Draw("nostack");
		for (auto &line : ex2_lines) line.Draw();
		c1->Print(pdf_file_name_2_q);
	}
	for (size_t i = 1; i < hist_ex2_q.size(); ++i) {
		hist_ex2_q[i].Draw();
        c1->Print(pdf_file_name_2_q);
	}
	c1->Print(pdf_file_name_2_q+"]");

	// close files
	opf.Close();
	ipf.Close();
	return 0;
}