#include <iostream>
#include <cmath>

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TString.h>

#include "include/event/ta_event.h"
#include "include/event/particle_type_event.h"

using namespace ribll;

constexpr double strip_pars[12][16][2] = {
	{
		{278, -0.00170},
		{279, -0.00169},
		{265, -0.00177},
		{264, -0.00172},
		{253, -0.00178},
		{248, -0.00184},
		{280, -0.00169},
		{240, -0.00181},
		{223, -0.00192},
		{232, -0.00189},
		{244, -0.00182},
		{199, -0.00204},
		{216, -0.00197},
		{205, -0.00227},
		{365, -0.00187},
		{319, -0.00163}
	},
	{
		{258, -0.00160},
		{253, -0.00175},
		{254, -0.00172},
		{252, -0.00164},
		{249, -0.00164},
		{234, -0.00180},
		{251, -0.00170},
		{248, -0.00170},
		{222, -0.00183},
		{238, -0.00187},
		{225, -0.00190},
		{197, -0.00208},
		{203, -0.00207},
		{259, -0.00208},
		{316, -0.00160},
		{308, -0.00161}
	},
	{
		{223, -0.00157},
		{251, -0.00164},
		{290, -0.00148},
		{260, -0.00162},
		{257, -0.00152},
		{251, -0.00175},
		{266, -0.00157},
		{258, -0.00162},
		{228, -0.00181},
		{229, -0.00180},
		{230, -0.00180},
		{212, -0.00185},
		{212, -0.00191},
		{292, -0.00159},
		{292, -0.00159},
		{292, -0.00159}
	},
	{
		{248, -0.00142},
		{279, -0.00153},
		{298, -0.00138},
		{307, -0.00137},
		{207, -0.00147},
		{265, -0.00156},
		{285, -0.00138},
		{293, -0.00142},
		{267, -0.00156},
		{267, -0.00156},
		{230, -0.00180},
		{247, -0.00177},
		{227, -0.00186},
		{193, -0.00229},
		{179, -0.00201},
		{333, -0.00141}
	},
	{
		{330, -0.00136},
		{328, -0.00132},
		{314, -0.00146},
		{331, -0.00134},
		{327, -0.00139},
		{310, -0.00147},
		{307, -0.00149},
		{284, -0.00155},
		{305, -0.00150},
		{296, -0.00151},
		{274, -0.00171},
		{278, -0.00174},
		{286, -0.00166},
		{264, -0.00175},
		{393, -0.00132},
		{393, -0.00132}
	},
	{
		{323, -0.00145},
		{341, -0.00135},
		{324, -0.00140},
		{311, -0.00143},
		{311, -0.00138},
		{259, -0.00152},
		{268, -0.00155},
		{286, -0.00139},
		{297, -0.00140},
		{246, -0.00149},
		{246, -0.00148},
		{256, -0.00158},
		{257, -0.00185},
		{257, -0.00185},
		{351, -0.00136},
		{351, -0.00136}
	},
	{
		{325, -0.00149},
		{330, -0.00138},
		{327, -0.00148},
		{320, -0.00149},
		{313, -0.00150},
		{298, -0.00160},
		{290, -0.00155},
		{279, -0.00147},
		{274, -0.00160},
		{276, -0.00163},
		{274, -0.00160},
		{249, -0.00171},
		{253, -0.00165},
		{242, -0.00180},
		{376, -0.00139},
		{375, -0.00139}
	},
	{
		{307, -0.00146},
		{291, -0.00144},
		{285, -0.00161},
		{295, -0.00154},
		{291, -0.00159},
		{282, -0.00154},
		{278, -0.00154},
		{276, -0.00157},
		{265, -0.00169},
		{241, -0.00177},
		{237, -0.00186},
		{227, -0.00203},
		{233, -0.00165},
		{357, -0.00147},
		{357, -0.00147}
	},
	{
		{363, -0.00119},
		{352, -0.00120},
		{357, -0.00132},
		{350, -0.00130},
		{352, -0.00130},
		{339, -0.00129},
		{331, -0.00142},
		{317, -0.00148},
		{272, -0.00181},
		{304, -0.00154},
		{293, -0.00154},
		{304, -0.00163},
		{280, -0.00170},
		{277, -0.00171},
		{427, -0.00123},
		{427, -0.00123},
	},
	{
		{310, -0.00129},
		{293, -0.00122},
		{327, -0.00128},
		{344, -0.00127},
		{349, -0.00127},
		{319, -0.00126},
		{311, -0.00137},
		{306, -0.00139},
		{300, -0.00143},
		{290, -0.00143},
		{283, -0.00147},
		{253, -0.00166},
		{260, -0.00174},
		{186, -0.00210},
		{391, -0.00122},
		{391, -0.00122}
	},
	{
		{439, -0.00122},
		{439, -0.00122},
		{439, -0.00122},
		{439, -0.00122},
		{372, -0.00119},
		{369, -0.00124},
		{369, -0.00129},
		{345, -0.00140},
		{313, -0.00143},
		{313, -0.00150},
		{309, -0.00149},
		{273, -0.00163},
		{264, -0.00170},
		{232, -0.00174},
		{439, -0.00122},
		{439, -0.00122}
	},
	{
		{435, -0.00129},
		{435, -0.00129},
		{435, -0.00129},
		{435, -0.00129},
		{373, -0.00119},
		{362, -0.00136},
		{338, -0.00151},
		{330, -0.00155},
		{322, -0.00145},
		{321, -0.00151},
		{316, -0.00152},
		{301, -0.00153},
		{297, -0.00155},
		{256, -0.00181},
		{435, -0.00129},
		{435, -0.00129}
	}
};


constexpr double pars[12][2] = {
	{319, -0.00163},
	{308, -0.00161},
	{292, -0.00159},
	{333, -0.00141},
	{393, -0.00132},
	{351, -0.00136},
	{376, -0.00139},
	{357, -0.00147},
	{427, -0.00123},
	{391, -0.00122},
	{439, -0.00122},
	{435, -0.00129}
};

double PidFit(double *x, double *par) {
	return (0.5/par[0]) * (
		sqrt(pow(x[0], 2.0) + 4.0*par[0]*pow(par[1]*x[0]-par[2], 2.0)) - x[0]
	);
}


int FitPid(int start_run, int end_run, int taf_index) {
	// input TAF telescope chain
	TChain chain("taf", "taf telescope");
	for (int run = start_run; run <= end_run; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		if (run > 716 && run < 739) continue;
		chain.AddFile(TString::Format(
			"%s%staf%d-telescope-ta-%04d.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			taf_index,
			run
		));
	}
	// input event
	TaEvent taf_event;
	// setup input branches
	taf_event.SetupInput(&chain);


	// input TAF pid chain
	TChain pid_chain("pid", "taf pid");
	for (int run = start_run; run <= end_run; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		if (run > 716 && run < 739) continue;
		pid_chain.AddFile(TString::Format(
			"%s%staf%d-particle-type-ta-%04d.root/tree",
			kGenerateDataPath,
			kParticleIdentifyDir,
			taf_index,
			run
		));
	}
	// add friend
	chain.AddFriend(&pid_chain);
	// input pid event
	ParticleTypeEvent pid_event;
	// setup input branches
	pid_event.SetupInput(&chain, "pid.");

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%staf%d-straight-pid-fit-%04d-%04d.root",
		kGenerateDataPath,
		kShowDir,
		taf_index,
		start_run,
		end_run
	);
	// output file
	TFile opf(output_file_name, "recreate");

	// graph for fitting
	// TAFD-CsI graph
	// 0: 1H, 1: 2H, 2: 3H
	TGraph taf_pid[2][3];
	TGraph taf_strip_pid[2][16][3];


	// total number of entries
	long long entries = chain.GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling TAF%d PID   0%%", taf_index);
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

		if (taf_event.num != 1) continue;
		if (pid_event.num != 1) continue;
		if (pid_event.charge[0] != 1) continue;

		// check CsI and get CsI index
		int csi_index = -1;
		if (taf_event.flag[0] == 0x3) csi_index = 0;
		else if (taf_event.flag[0] == 0x5) csi_index = 1;
		if (csi_index < 0) continue;

		// fill to graph
		if (pid_event.mass[0] >= 1 && pid_event.mass[0] <= 3) {
			int mass_index = pid_event.mass[0]-1;
			taf_pid[csi_index][mass_index].AddPoint(
				taf_event.energy[0][1],
				taf_event.energy[0][0] * cos(taf_event.theta[0])
			);
			int fs = taf_event.front_strip[0];
			taf_strip_pid[csi_index][fs][mass_index].AddPoint(
				taf_event.energy[0][1],
				taf_event.energy[0][0]
			);
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fit TAFD-CsI PID graph
	for (int i = 0; i < 2; ++i) {
		std::cout << "TAFD-CsI" << i << " straight PID parameters.\n";
		for (int j = 0; j < 3; ++j) {
			TF1 *f1 = new TF1(
				TString::Format("fc%dh%d", i, j+1), PidFit, 0.0, 2e4, 3
			);
			f1->SetParameter(0, 250.0);
			f1->SetParameter(1, -0.0014);
			f1->SetParameter(2, 50.0);
			f1->SetParLimits(2, 0.0, 1e10);
			taf_pid[i][j].Fit(f1, "RQ+");
			std::cout << f1->GetParameter(0)
				<< ", " << f1->GetParameter(1)
				<< ", " << f1->GetParameter(2)
				<< "\n";
		}
	}

	// fit TAFD-CsI strip PID graph
	for (int i = 0; i < 2; ++i) {
		std::cout << "TAFD-CsI" << i << " strip straight PID parameters.\n";
		for (int j = 0; j < 16; ++j) {
			std::cout << "Strip " << j << "\n";
			for (int k = 0; k < 3; ++k) {
				TF1 *f1 = new TF1(
					TString::Format("fc%ds%dh%d", i, j, k+1), PidFit, 0.0, 2e4, 3
				);
				f1->SetParameter(0, 300.0);
				f1->SetParameter(1, -0.0015);
				f1->SetParameter(2, 80.0);
				f1->SetParLimits(2, 0.0, 1e10);
				taf_strip_pid[i][j][k].Fit(f1, "RQ+");
				std::cout << f1->GetParameter(0)
					<< ", " << f1->GetParameter(1)
					<< ", " << f1->GetParameter(2)
					<< "\n";
			}
		}
	}


	// save pid graph
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 3; ++j) {
			taf_pid[i][j].Write(TString::Format(
				"gc%dh%d", i, j+1
			));
		}
	}
	// save pid strip graph
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 16; ++j) {
			for (int k = 0; k < 3; ++k) {
				taf_strip_pid[i][j][k].Write(TString::Format(
					"gc%ds%dh%d", i, j, k+1
				));
			}
		}
	}
	// close files
	opf.Close();

	return 0;
}


double Gaus3(double *x, double *par) {
	return
		par[0]*exp(-0.5*pow((x[0]-par[1])/par[2], 2.0))
		+ par[3]*exp(-0.5*pow((x[0]-par[4])/par[5], 2.0))
		+ par[6]*exp(-0.5*pow((x[0]-par[7])/par[8], 2.0));
}

int StraightPid(int start_run, int end_run) {
	// input taf chain
	TChain chain("taf0", "taf telescope");
	// input event
	TaEvent taf_event[6];
	// add file
	for (int run = start_run; run <= end_run; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		if (run > 716 && run < 739) continue;
		chain.AddFile(TString::Format(
			"%s%staf0-telescope-ta-%04d.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			run
		));
	}
	// setup input branches
	taf_event[0].SetupInput(&chain);	

	// add TAF1-TAF5
	for (int taf_index = 1; taf_index < 6; ++taf_index) {
		TChain *taf_chain = new TChain(
			TString::Format("taf%d", taf_index), "taf telescope"
		);
		for (int run = start_run; run <= end_run; ++run) {
			if (run == 628) continue;
			if (run > 652 && run < 675) continue;
			if (run > 716 && run < 739) continue;
			taf_chain->AddFile(TString::Format(
				"%s%staf%d-telescope-ta-%04d.root/tree",
				kGenerateDataPath,
				kTelescopeDir,
				taf_index,
				run
			));
		}
		// add friend
		chain.AddFriend(taf_chain);
		// setup input branches
		taf_event[taf_index].SetupInput(
			&chain, TString::Format("taf%d.", taf_index).Data()
		);
	}

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%staf-straight-pid-%04d-%04d.root",
		kGenerateDataPath,
		kShowDir,
		start_run,
		end_run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// straight PID
	TH2F hist_straight_pid[12];
	for (int i = 0; i < 12; ++i) {
		hist_straight_pid[i] = TH2F(
			TString::Format("hsc%d", i),
			"straight PID",
			2000, 0, 20000, 1200, 0, 120 
		);
	}
	// projected corrected energy
	TH1F hist_ef[12];
	for (int i = 0; i < 12; ++i) {
		hist_ef[i] = TH1F(
			TString::Format("hefc%d", i),
			"particle fixed energy",
			1200, 0, 120
		);
	}
	// straight PID for single strip
	TH2F hist_strip_straight_pid[12][16];
	for (int i = 0; i < 12; ++i) {
		for (int j = 0; j < 16; ++j) {
			hist_strip_straight_pid[i][j] = TH2F(
				TString::Format("hsc%ds%d", i, j),
				"straight PID for single strip",
				2000, 0, 20000, 300, 0, 120
			);
		}
	}
	// projected fixed energy for single strip
	TH1F hist_strip_ef[12][16];
	for (int i = 0; i < 12; ++i) {
		for (int j = 0; j < 16; ++j) {
			hist_strip_ef[i][j] = TH1F(
				TString::Format("hefc%ds%d", i, j),
				"fixed energy for single strip",
				300, 0, 120
			);
		}
	}


	// total number of entries
	long long entries = chain.GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling TAF straight PID   0%%");
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

		for (int taf_index = 0; taf_index < 6; ++taf_index) {
			if (taf_event[taf_index].num != 1) continue;

			// check CsI and get CsI index
			int csi_index = -1;
			if (taf_event[taf_index].flag[0] == 0x3) csi_index = 0;
			else if (taf_event[taf_index].flag[0] == 0x5) csi_index = 1;
			if (csi_index < 0) continue;
			int fs = taf_event[taf_index].front_strip[0];

			// whole CsI
			double cde =
				taf_event[taf_index].energy[0][0]
				* cos(taf_event[taf_index].theta[0]);
			double e = taf_event[taf_index].energy[0][1];
			double a = pars[taf_index*2+csi_index][0];
			double b = pars[taf_index*2+csi_index][1];
			double ef = sqrt(cde*e + a*cde*cde) + b*e;
			hist_straight_pid[taf_index*2+csi_index].Fill(e, ef);
			hist_ef[taf_index*2+csi_index].Fill(ef);
			
			// single strip
			double de = taf_event[taf_index].energy[0][0];
			double as = strip_pars[taf_index*2+csi_index][fs][0];
			double bs = strip_pars[taf_index*2+csi_index][fs][1];
			double efs = sqrt(de*e + as*de*de) + bs*e;
			hist_strip_straight_pid[taf_index*2+csi_index][fs].Fill(e, efs);
			hist_strip_ef[taf_index*2+csi_index][fs].Fill(efs);

		}

	}
	// show finish
	printf("\b\b\b\b100%%\n");

	double init_pars[9] = {
		500, 55, 3,
		500, 75, 3,
		200, 85, 3
	};
	for (int i = 0; i < 12; ++i) {
		// fit 1H, 2H and 3H
		TF1 *f1 = new TF1(TString::Format("fc%d", i), Gaus3, 30, 120, 9);
		f1->SetParameters(init_pars);
		f1->SetParLimits(1, 40, 60);
		f1->SetParLimits(2, 1, 10);
		f1->SetParLimits(4, 60, 80);
		f1->SetParLimits(5, 1, 10);
		f1->SetParLimits(7, 70, 90);
		f1->SetParLimits(8, 1, 10);
		if (i < 3) {
			f1->SetParameter(1, 50);
			f1->SetParameter(4, 70);
			f1->SetParameter(7, 80);
		}
		hist_ef[i].Fit(f1, "RQ+");
		for (int j = 0; j < 3; ++j) {
			TF1 *f2 = new TF1(
				TString::Format("f%dc%d", j, i), "gaus", 30, 120
			);
			f2->SetLineColor(kGreen);
			f2->SetParameters(f1->GetParameters()+j*3);
			hist_ef[i].GetListOfFunctions()->Add(f2);
		}
		// std::cout << "TAF-CsI " << i << ":\n"
		// 	<< "  " << f1->GetParameter(1)
		// 	<< ", " << f1->GetParameter(2) << "\n"
		// 	<< "  " << f1->GetParameter(4)
		// 	<< ", " << f1->GetParameter(5) << "\n"
		// 	<< "  " << f1->GetParameter(7)
		// 	<< ", " << f1->GetParameter(8) << "\n"
		// 	<< "  " << (f1->GetParameter(7)-f1->GetParameter(4))
		// 		/ (f1->GetParameter(8)+f1->GetParameter(5))
		// 	<< "\n";
		std::cout << f1->GetParameter(4) << ", "
			<< f1->GetParameter(5) << "\n";
	}

	for (int i = 0; i < 12; ++i) {
		std::cout << "TAF-CsI " << i << " strips\n";
		for (int j = 0; j < 16; ++j) {
			// fit 1H, 2H and 3H
			TF1 *f1 = new TF1(
				TString::Format("fc%ds%d", i, j), Gaus3, 30, 120, 9
			);
			f1->SetParameters(init_pars);
			f1->SetParameter(0, 200);
			f1->SetParameter(3, 100);
			f1->SetParameter(6, 100);
			f1->SetParLimits(1, 40, 70);
			f1->SetParLimits(2, 1, 10);
			f1->SetParLimits(4, 60, 95);
			f1->SetParLimits(5, 1, 10);
			f1->SetParLimits(7, 70, 110);
			f1->SetParLimits(8, 1, 10);
			if (i < 2) {
				f1->SetParameter(1, 50);
				f1->SetParameter(4, 70);
				f1->SetParameter(7, 80);
				if (j > 5) {
					f1->SetParameter(1, 60);
					f1->SetParameter(4, 75);
					f1->SetParameter(7, 90);
				}
			} else if (i == 2) {
				if (j==6 || j==7 || j==8 || j==10 || j==11) {
					f1->SetParameter(1, 60);
					f1->SetParameter(4, 80);
					f1->SetParameter(7, 90);
				} else if (j == 13) {
					f1->SetParameter(1, 65);
					f1->SetParameter(4, 85);
					f1->SetParameter(7, 95);
				}
			} else if (i == 3) {
				if (j!=0 && j!=1 && j!=4 && j!=8 && j!=9 && j!=13) {
					f1->SetParameter(1, 65);
					f1->SetParameter(4, 85);
					f1->SetParameter(7, 95);
				} else if (j == 13) {
					f1->SetParameter(1, 60);
					f1->SetParameter(4, 80);
					f1->SetParameter(7, 90);
				}
			} else if (i >= 4 && i <= 9) {
				f1->SetParameter(1, 60);
				f1->SetParameter(4, 80);
				f1->SetParameter(7, 95);
				if (!(i == 7 && j == 10) && j >= 8) {
					f1->SetParameter(1, 65);
					f1->SetParameter(4, 85);
					f1->SetParameter(7, 100);
				}
			} else if (i == 10 || i == 11) {
				f1->SetParameter(1, 65);
				f1->SetParameter(4, 85);
				f1->SetParameter(7, 100);
			}

			hist_strip_ef[i][j].Fit(f1, "RQ+");
			for (int k = 0; k < 3; ++k) {
				TF1 *f2 = new TF1(
					TString::Format("f%dc%ds%d", j, i, k), "gaus", 30, 120
				);
				f2->SetLineColor(kGreen);
				f2->SetParameters(f1->GetParameters()+k*3);
				hist_strip_ef[i][j].GetListOfFunctions()->Add(f2);
			}
			std::cout << f1->GetParameter(1) << ", "
				<< f1->GetParameter(2) << "\n";
		}
	}

	// save histograms
	for (int i = 0; i < 12; ++i) {
		hist_straight_pid[i].Write();
		hist_ef[i].Write();
	}
	for (int i = 0; i < 12; ++i) {
		for (int j = 0; j < 16; ++j) {
			hist_strip_straight_pid[i][j].Write();
			hist_strip_ef[i][j].Write();
		}
	}
	// close files
	opf.Close();

	return 0;
}


/// @brief print program usage
/// @param[in] name program name 
void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] start_run end_run taf_index\n"
		"  start_run         Set run number.\n"
		"  end_run           Set the last run to chain, included.\n"
		"  taf_index         Set TAF index, 0-5.\n"
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
	if (pos_start+1 >= argc || (fit && pos_start+2 >= argc)) {
		// positional arguments less than 3
		std::cerr << "Error: Miss detector argument.\n";
		PrintUsage(argv[0]);
		return -1;
	}
	// start run number
	int start_run = atoi(argv[pos_start]);
	// end run number
	int end_run = atoi(argv[pos_start+1]);
	// taf index
	int taf_index = fit ? atoi(argv[pos_start+2]) : -1;

	if (fit) {
		if (FitPid(start_run, end_run, taf_index)) {
			return -1;
		}
	} else {
		if (StraightPid(start_run, end_run)) {
			return -1;
		}
	}

	return 0;
}