#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TString.h>

#include "include/event/t0_event.h"

using namespace ribll;

double PidFit(double *x, double *par) {
	return (0.5/par[0]) * (
		sqrt(pow(x[0], 2.0) + 4.0*par[0]*pow(par[1]*x[0]-par[2], 2.0)) - x[0]
	);
}


int FitPid(int start_run, int end_run) {
	// input t0 chain
	TChain chain("t0", "t0 telescope");
	for (int run = start_run; run <= end_run; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		if (run > 716 && run < 739) continue;
		chain.AddFile(TString::Format(
			"%s%st0-telescope-ta-v3-%04d.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			run
		));
	}
	// input event
	T0Event t0_event;
	// setup input branches
	t0_event.SetupInput(&chain);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%st0-straight-pid-fit-%04d-%04d.root",
		kGenerateDataPath,
		kShowDir,
		start_run,
		end_run
	);
	// output file
	TFile opf(output_file_name, "recreate");

	// graph for fitting
	// D1D2 PID graph
	// 0: 4He, 1: 6Li, 2: 7Li, 3: 9Be, 4: 10Be, 5: 11B, 6: 12B
	TGraph d1d2_pid[7];
	// D2D3 PID graph
	// 0: 4He, 1: 6Li, 2: 7Li, 3: 9Be, 4: 10Be
	TGraph d2d3_pid[5];


	// total number of entries
	long long entries = chain.GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling T0 PID   0%%");
	fflush(stdout);
	// loop to fill pid histogram
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		chain.GetEntry(entry);

		for (int i = 0; i < t0_event.num; ++i) {
			if (t0_event.hole[i]) continue;
			if (t0_event.layer[i] == 1) {
				int index = -1;
				if (t0_event.charge[i] == 2 && t0_event.mass[i] == 4) {
					index = 0;
				} else if (t0_event.charge[i] == 3 && t0_event.mass[i] == 6) {
					index = 1;
				} else if (t0_event.charge[i] == 3 && t0_event.mass[i] == 7) {
					index = 2;
				} else if (t0_event.charge[i] == 4 && t0_event.mass[i] == 9) {
					index = 3;
				} else if (t0_event.charge[i] == 4 && t0_event.mass[i] == 10) {
					index = 4;
				} else if (t0_event.charge[i] == 5 && t0_event.mass[i] == 11) {
					index = 5;
				} else if (t0_event.charge[i] == 5 && t0_event.mass[i] == 12) {
					index = 6;
				}
				if (index >= 0 && index <= 6) {
					d1d2_pid[index].AddPoint(
						t0_event.energy[i][1], t0_event.energy[i][0]
					);
				}
			} else if (t0_event.layer[i] == 2) {
				int index = -1;
				if (t0_event.charge[i] == 2 && t0_event.mass[i] == 4) {
					index = 0;
				} else if (t0_event.charge[i] == 3 && t0_event.mass[i] == 6) {
					index = 1;
				} else if (t0_event.charge[i] == 3 && t0_event.mass[i] == 7) {
					index = 2;
				} else if (t0_event.charge[i] == 4 && t0_event.mass[i] == 9) {
					index = 3;
				} else if (t0_event.charge[i] == 4 && t0_event.mass[i] == 10) {
					index = 4;
				}
				if (index >= 0 && index <= 4) {
					d2d3_pid[index].AddPoint(
						t0_event.energy[i][2], t0_event.energy[i][1]
					);
				}
			}
		}		
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fit D1D2 PID graph
	// D1D2 fit range
	double d1d2_range[7] = {1.2e4, 2e4, 2.2e4, 3.2e4, 3.5e4, 4.2e4, 4.5e4};
	// D1D2 initialize energy-fixed-paramter (c)
	double d1d2_init_c[7] = {6.4e3, 1.2e4, 1.2e4, 2e4, 2e4, 2.8e4, 2.8e4};
	std::cout << "D1D2 straight PID parameters.\n";
	for (int i = 0; i < 7; ++i) {
		TF1 *f1 = new TF1(
			TString::Format("fd1d2%d", i), PidFit, 0.0, d1d2_range[i], 3
		);
		f1->SetParameter(0, 0.47);
		f1->SetParameter(1, -0.042);
		f1->SetParameter(2, d1d2_init_c[i]);
		f1->SetParLimits(2, 0.0, 1e10);
		d1d2_pid[i].Fit(f1, "RQ+");
		std::cout << "A " << f1->GetParameter(0)
			<< ", B " << f1->GetParameter(1)
			<< ", C " << f1->GetParameter(2) << "\n";
	}

	// fit D1D2 PID graph
	// D1D2 fit range
	double d2d3_range[5] = {1.2e4, 2.2e4, 2.3e4, 2.6e4, 3e4};
	// D1D2 initialize energy-fixed-paramter (c)
	double d2d3_init_c[5] = {8e3, 1.3e4, 1.5e4, 2.5e4, 2.5e4};
	std::cout << "D2D3 straight PID parameters.\n";
	for (int i = 0; i < 5; ++i) {
		TF1 *f1 = new TF1(
			TString::Format("fd2d3%d", i), PidFit, 0.0, d2d3_range[i], 3
		);
		f1->SetParameter(0, 0.54);
		f1->SetParameter(1, -0.042);
		f1->SetParameter(2, d2d3_init_c[i]);
		f1->SetParLimits(2, 0.0, 1e10);
		d2d3_pid[i].Fit(f1, "RQ+");
		std::cout << "A " << f1->GetParameter(0)
			<< ", B " << f1->GetParameter(1)
			<< ", C " << f1->GetParameter(2) << "\n";
	}

	// save pid histograms
	d1d2_pid[0].Write("gd1d2he4");
	d1d2_pid[1].Write("gd1d2li6");
	d1d2_pid[2].Write("gd1d2li7");
	d1d2_pid[3].Write("gd1d2be9");
	d1d2_pid[4].Write("gd1d2be10");
	d1d2_pid[5].Write("gd1d2b11");
	d1d2_pid[6].Write("gd1d2b12");
	d2d3_pid[0].Write("gd2d3he4");
	d2d3_pid[1].Write("gd2d3li6");
	d2d3_pid[2].Write("gd2d3li7");
	d2d3_pid[3].Write("gd2d3be9");
	d2d3_pid[4].Write("gd2d3be10");
	// close files
	opf.Close();

	return 0;
}


int StraightPid(int start_run, int end_run, double *parameters) {
	// input t0 chain
	TChain chain("t0", "t0 telescope");
	for (int run = start_run; run <= end_run; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		if (run > 716 && run < 739) continue;
		chain.AddFile(TString::Format(
			"%s%st0-telescope-ta-v1-%04d.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			run
		));
	}
	// input event
	T0Event t0_event;
	// setup input branches
	t0_event.SetupInput(&chain);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%st0-straight-pid-%04d-%04d.root",
		kGenerateDataPath,
		kShowDir,
		start_run,
		end_run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// D1D2 PID
	TH2F d1d2_pid("hpd1d2", "D1D2 PID", 2500, 0, 50000, 3000, 0, 60000);
	// D2D3 PID
	TH2F d2d3_pid("hpd2d3", "D2D3 PID", 1750, 0, 35000, 2000, 0, 40000);
	// straight D1D2 PID
	TH2F d1d2_straight_pid(
		"hspd1d2", "D1D2 straight PID", 2500, 0, 50000, 3000, 0, 60000
	);
	// straight D2D3 PID
	TH2F d2d3_straight_pid(
		"hspd2d3", "D2D3 straight PID", 1750, 0, 35000, 2000, 0, 40000
	);
	// projected D1D2 corrected energy
	TH1F d1d2ef("hd1d2ef", "D1D2 particle fixed energy", 3000, 0, 60000);
	// projected D2D3 corrected energy
	TH1F d2d3ef("hd2d3ef", "D2D3 particle fixed energy", 2000, 0, 40000); 

	double a12 = parameters[0];
	double b12 = parameters[1];
	double a23 = parameters[2];
	double b23 = parameters[3];

	// total number of entries
	long long entries = chain.GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling T0 straight PID   0%%");
	fflush(stdout);
	// loop to fill pid histogram
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		chain.GetEntry(entry);

		for (int i = 0; i < t0_event.num; ++i) {
			if (t0_event.hole[i]) continue;
			if (t0_event.flag[i] == 0x3) {
				double de = t0_event.energy[i][0];
				double e = t0_event.energy[i][1];
				double ef = sqrt(de*e + a12*de*de) + b12*e;
				d1d2_pid.Fill(e, de);
				d1d2_straight_pid.Fill(e, ef);
				if (e > 1000) d1d2ef.Fill(ef);
			} else if (t0_event.flag[i] == 0x7 && t0_event.ssd_flag == 0) {
				double de = t0_event.energy[i][1];
				double e = t0_event.energy[i][2];
				double ef = sqrt(de*e + a23*de*de) + b23*e;
				d2d3_pid.Fill(e, de);
				d2d3_straight_pid.Fill(e, ef);
				if (e > 500) d2d3ef.Fill(ef);
			}
		}		
	}
	// show finish
	printf("\b\b\b\b100%%\n");
	
	// save histograms
	d1d2_pid.Write();
	d2d3_pid.Write();
	d1d2_straight_pid.Write();
	d2d3_straight_pid.Write();
	d1d2ef.Write();
	d2d3ef.Write();
	// close files
	opf.Close();

	return 0;
}


/// @brief print program usage
/// @param[in] name program name 
void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] start_run end_run\n"
		"  start_run         Set run number.\n"
		"  end_run           Set the last run to chain, included.\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -f                Fit PID and get straight parameters.\n";
}

/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] fit fit and get parameters
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	bool &fit
) {
	// initialize
	help = false;
	fit = false;
	// start index of positional arugments
	int result = 0;
	for (result = 1; result < argc; ++result) {
		// assumed that all options have read
		if (argv[result][0] != '-') break;
		// short option contains only one letter
		if (argv[result][2] != 0) return -result;
		if (argv[result][1] == 'h') {
			help = true;
			return result;
		} else if (argv[result][1] == 'f') {
			fit = true;
		} else {
			return -result;
		}
	}
	return result;
}


int main(int argc, char **argv) {
	if (argc < 3) {
		PrintUsage(argv[0]);
		return -1;
	}
	bool help = false;
	bool fit = false;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(
		argc, argv, help, fit
	);
	// need help
	if (help) {
		PrintUsage(argv[0]);
		return 0;
	}
	// invalid arguments
	if (pos_start < 0) {
		if (-pos_start < argc) {
			std::cerr << "Error: Invaild option " << argv[-pos_start] << ".\n";
		} else {
			std::cerr << "Error: Option need parameter.\n";
		}
		PrintUsage(argv[0]);
		return -1;
	}
	// check number of positional arguments
	if (pos_start+1 >= argc) {
		// positional arguments less than 3
		std::cerr << "Error: Miss detector argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}
	// start run number
	int start_run = atoi(argv[pos_start]);
	// end run number
	int end_run = atoi(argv[pos_start+1]);

	double parameters[] = {
		0.47, -0.042,
		0.59, -0.042
	};

	if (fit) {
		if (FitPid(start_run, end_run)) {
			return -1;
		}
	} else {
		if (StraightPid(start_run, end_run, parameters)) {
			return -1;
		}
	}

	return 0;
}