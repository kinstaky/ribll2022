#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

#include "include/defs.h"

using namespace ribll;

constexpr int bins = 50;
constexpr double step = 0.1;

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
	TH1F hist_ex0_qs[5][8][8];
	TH1F hist_ex1_qs[5][8][8];
	TH1F hist_ex2_qs[5][8][8];
	// ground state excited energy spectrum, fix Q, change start point
	TH1F hist_ex0_s[4][8];
	TH1F hist_ex1_s[4][8];
	TH1F hist_ex2_s[4][8];

	// set histograms parameters, changes Q range and start point
	// ground state
	double ex0_qs_range[5];
	double ex0_qs_min[5][8];
	double ex0_qs_max[5][8];
	for (int i = 0; i < 5; ++i) {
		ex0_qs_range[i] = 4.0 - i * 0.5;
		for (int j = 0; j < i+1; ++j) {
			ex0_qs_min[i][j] = -14.0 + 0.5*j;
			ex0_qs_max[i][j] = ex0_qs_min[i][j] + ex0_qs_range[i];
		}
	}
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < i+1; ++j) {
			for (int k = 0; k < 6; ++k) {
				hist_ex0_qs[i][j][k] = TH1F(
					TString::Format("hex0r%dq%de%d", i, j, k),
					TString::Format(
						"ex0 Q %lf to %lf, ex %lf to %lf",
						ex0_qs_min[i][j],
						ex0_qs_min[i][j]+ex0_qs_range[i],
						12.0+0.05*k,
						27.0+0.05*k
					),
					bins, 12.0+0.05*k, 27.0+0.05*k
				);
			}
		}
	}
	// first excited state
	double ex1_qs_range[5];
	double ex1_qs_min[5][8];
	double ex1_qs_max[5][8];
	for (int i = 0; i < 5; ++i) {
		ex1_qs_range[i] = 3.0 - i * 0.5;
		for (int j = 0; j < i+1; ++j) {
			ex1_qs_min[i][j] = -17.0 + 0.5*j;
			ex1_qs_max[i][j] = ex1_qs_min[i][j] + ex1_qs_range[i];
		}
	}
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < i+1; ++j) {
			for (int k = 0; k < 6; ++k) {
				hist_ex1_qs[i][j][k] = TH1F(
					TString::Format("hex1r%dq%de%d", i, j, k),
					TString::Format(
						"ex1 Q %lf to %lf, ex %lf to %lf",
						ex1_qs_min[i][j],
						ex1_qs_min[i][j]+ex1_qs_range[i],
						12.0+0.05*k,
						27.0+0.05*k
					),
					bins, 12.0+0.05*k, 27.0+0.05*k
				);
			}
		}
	}
	// 6MeV state
	double ex2_qs_range[5];
	double ex2_qs_min[5][8];
	double ex2_qs_max[5][8];
	for (int i = 0; i < 5; ++i) {
		ex2_qs_range[i] = 3.0 - i * 0.5;
		for (int j = 0; j < i+1; ++j) {
			ex2_qs_min[i][j] = -20.0 + 0.5*j;
			ex2_qs_max[i][j] = ex2_qs_min[i][j] + ex2_qs_range[i];
		}
	}
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < i+1; ++j) {
			for (int k = 0; k < 6; ++k) {
				hist_ex2_qs[i][j][k] = TH1F(
					TString::Format("hex2r%dq%de%d", i, j, k),
					TString::Format(
						"ex2 Q %lf to %lf, ex %lf to %lf",
						ex2_qs_min[i][j],
						ex2_qs_min[i][j]+ex2_qs_range[i],
						12.0+0.05*k,
						27.0+0.05*k
					),
					bins, 12.0+0.05*k, 27.0+0.05*k
				);
			}
		}
	}

	// set histograms parameters, fix Q range and change start point
	// ground state
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 8; ++j) {
			hist_ex0_s[i][j] = TH1F(
				TString::Format("hex0b%de%d", i, j),
				TString::Format(
					"ex0 Q -13.0 to -11.0, %d bins, ex %lf to %lf",
					50 + 10 * i,
					12.0 + 0.05 * j,
					27.0 + 0.05 * j
				),
				50+10*i, 12.0+0.05*j, 27.0+0.05*j
			);
		}
	}
	// // first excited state
	// for (int i = 0; i < 4; ++i) {
	// 	for (int j = 0; j < 8; ++j) {
	// 		hist_ex1_s[i][j] = TH1F(
	// 			TString::Format("hex1b%de%d", i, j),
	// 			TString::Format(
	// 				"ex1 Q -16.5 to -14.5, %d bins, ex %lf to %lf",
	// 				50 + 10 * i,
	// 				12.0 + 0.05 * j,
	// 				27.0 + 0.05 * j
	// 			),
	// 			50+10*i, 12.0+0.05*j, 27.0+0.05*j
	// 		);
	// 	}
	// }
	// // 6MeV state
	// for (int i = 0; i < 4; ++i) {
	// 	for (int j = 0; j < 8; ++j) {
	// 		hist_ex2_s[i][j] = TH1F(
	// 			TString::Format("hex2b%de%d", i, j),
	// 			TString::Format(
	// 				"ex0 Q -13.0 to -11.0, %d bins, ex %lf to %lf",
	// 				50 + 10 * i,
	// 				12.0 + 0.05 * j,
	// 				27.0 + 0.05 * j
	// 			),
	// 			50+10*i, 12.0+0.05*j, 27.0+0.05*j
	// 		);
	// 	}
	// }


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
		for (int i = 0; i < 5; ++i) {
			for (int j = 0; j < i+1; ++j) {
				if (q[0] < ex0_qs_min[i][j]) continue;
				if (q[0] > ex0_qs_max[i][j]) continue;
				for (int k = 0; k < 6; ++k) {
					hist_ex0_qs[i][j][k].Fill(stateless_excited_energy[0][0]);
				}
			}
		}
		// first excited state
		for (int i = 0; i < 5; ++i) {
			for (int j = 0; j < i+1; ++j) {
				if (q[0] < ex1_qs_min[i][j]) continue;
				if (q[0] > ex1_qs_max[i][j]) continue;
				for (int k = 0; k < 6; ++k) {
					hist_ex1_qs[i][j][k].Fill(stateless_excited_energy[1][0]);
				}
			}
		}
		// 6MeV state
		for (int i = 0; i < 5; ++i) {
			for (int j = 0; j < i+1; ++j) {
				if (q[0] < ex2_qs_min[i][j]) continue;
				if (q[0] > ex2_qs_max[i][j]) continue;
				for (int k = 0; k < 6; ++k) {
					hist_ex2_qs[i][j][k].Fill(stateless_excited_energy[2][0]);
				}
			}
		}

		// ground state
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 8; ++j) {
				if (q[0] < -13.0 || q[0] > -11.0) continue;
				hist_ex0_s[i][j].Fill(stateless_excited_energy[0][0]);
			}
		}
	}

	// save images
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < i+1; ++j) {
			for (int k = 0; k < 6; ++k) {
				hist_ex0_qs[i][j][k].Write();
			}
		}
	}
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < i+1; ++j) {
			for (int k = 0; k < 6; ++k) {
				hist_ex1_qs[i][j][k].Write();
			}
		}
	}
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < i+1; ++j) {
			for (int k = 0; k < 6; ++k) {
				hist_ex2_qs[i][j][k].Write();
			}
		}
	}
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 8; ++j) {
			hist_ex0_s[i][j].Write();
		}
	}

	// save as pdf file, ground state
	TCanvas *c1 = new TCanvas("c1", "c1", 1920, 1080);
	c1->cd();
	TString pdf_file_name_0 = TString::Format(
		"%simage/show-spectrum2%s-ex0-qs.pdf",
		kGenerateDataPath,
		suffix.empty() ? "" : ("-"+suffix).c_str()
	);
	c1->Print(pdf_file_name_0+"[");
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < i+1; ++j) {
			for (int k = 0; k < 6; ++k) {
				hist_ex0_qs[i][j][k].Draw();
				c1->Print(pdf_file_name_0);
			}
		}
	}
	c1->Print(pdf_file_name_0+"]");
	// first excited state
	c1->cd();
	TString pdf_file_name_1 = TString::Format(
		"%simage/show-spectrum2%s-ex1-qs.pdf",
		kGenerateDataPath,
		suffix.empty() ? "" : ("-"+suffix).c_str()
	);
	c1->Print(pdf_file_name_1+"[");
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < i+1; ++j) {
			for (int k = 0; k < 6; ++k) {
				hist_ex1_qs[i][j][k].Draw();
				c1->Print(pdf_file_name_1);
			}
		}
	}
	c1->Print(pdf_file_name_1+"]");
	// 6MeV state
	c1->cd();
	TString pdf_file_name_2 = TString::Format(
		"%simage/show-spectrum2%s-ex2-qs.pdf",
		kGenerateDataPath,
		suffix.empty() ? "" : ("-"+suffix).c_str()
	);
	c1->Print(pdf_file_name_2+"[");
	for (int i = 0; i < 5; ++i) {
		for (int j = 0; j < i+1; ++j) {
			for (int k = 0; k < 6; ++k) {
				hist_ex2_qs[i][j][k].Draw();
				c1->Print(pdf_file_name_2);
			}
		}
	}
	c1->Print(pdf_file_name_2+"]");

	// ground state, -13.0 to -11.0 MeV
	c1->cd();
	TString pdf_file_name_0s = TString::Format(
		"%simage/show-spectrum2%s-ex0-s.pdf",
		kGenerateDataPath,
		suffix.empty() ? "" : ("-"+suffix).c_str()
	);
	c1->Print(pdf_file_name_0s+"[");
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 8; ++j) {
			hist_ex0_s[i][j].Draw();
			c1->Print(pdf_file_name_0s);
		}
	}
	c1->Print(pdf_file_name_0s+"]");


	// close files
	opf.Close();
	ipf.Close();
	return 0;
}