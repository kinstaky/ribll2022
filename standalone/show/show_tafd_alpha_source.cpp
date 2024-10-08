#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

#include "include/event/dssd_event.h"

using namespace ribll;

// run used for alpha calibration
const unsigned int alpha_calibration_run[6] = {816, 816, 825, 825, 825, 825};

int main(int argc, char **argv) {
	if (argc != 2) {
		std::cout << "Usage: " << argv[0] << " taf_index\n";
		return -1;
	}

	int taf_index = atoi(argv[1]);

	// calibrate tafd with alpha source
	// output root file name
	TString output_file_name;
	output_file_name.Form(
		"%s%stafd%d-alpha-peaks.root",
		kGenerateDataPath,
		kShowDir,
		taf_index
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output energy histogram for each strip
	TH1F hist_energy[16];
	TH1F hist_channel[16];
	for (int i = 0; i < 16; ++i) {
		hist_channel[i] = TH1F(
			TString::Format("hc%d", i),
			TString::Format("before calibration, strip %d", i),
			500,
			taf_index < 2 ? 1500 : 11500,
			taf_index < 2 ? 2500 : 17500
		);
		hist_energy[i] = TH1F(
			TString::Format("he%d", i),
			TString::Format("after calibration, strip %d", i),
			taf_index < 2 ? 100 : 500,
			4.5, 6.5
		);
	}
	if (taf_index == 5) {
		hist_channel[15].SetTitle("bad strip");
		hist_energy[15].SetTitle("bad strip");
	}

	// parameters txt file
	TString param_file_name = TString::Format(
		"%s%stafd%d-alpha-cali-param-nodead.txt",
		kGenerateDataPath,
		kCalibrationDir,
		taf_index
	);
	std::ifstream fin(param_file_name.Data());
	if (!fin.good()) {
		std::cerr << "Error: Open file "
			<< param_file_name << " failed.\n";
		return -1;
	}
	double cali_p0[16], cali_p1[16], tmp;
	for (int i = 0; i < 16; ++i) {
		fin >> cali_p0[i] >> cali_p1[i] >> tmp >> tmp >> tmp;
	}
	fin.close();

	if (taf_index < 2) {
		// VME
		// input file name
		TString input_file_name;
		input_file_name.Form(
			"%s%s%04u.root",
			kCrate3Path,
			kCrate3FileName,
			alpha_calibration_run[taf_index]
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
		// energy values
		int madc[2][32];
		// setup input branches
		ipt->SetBranchAddress("madc", madc);

		// madc module of this tafd
		size_t mod = vtaf_front_module[taf_index];
		// madc channel of this tafd
		size_t ch = vtaf_front_channel[taf_index];

		// total number of entries
		long long entries = ipt->GetEntries();
		// 1/100 of entries
		long long entry100 = entries / 100 + 1;
		// show start
		printf("Filling histogram   0%%");
		fflush(stdout);
		// fill energy to histogram
		for (long long entry = 0; entry < entries; ++entry) {
			// show process
			if (entry % entry100 == 0) {
				printf("\b\b\b\b%3lld%%", entry / entry100);
				fflush(stdout);
			}
			ipt->GetEntry(entry);
			for (size_t i = 0; i < 16; ++i) {
				if (madc[mod][ch+i] > 1000 && madc[mod][ch+i] < 5000) {
					hist_channel[i].Fill(madc[mod][ch+i]);
					hist_energy[i].Fill(cali_p0[i] + cali_p1[i]*madc[mod][ch+i]);
				}
			}
		}
		// show finish
		printf("\b\b\b\b100%%\n");
		// close input file
		ipf.Close();

	} else {
		// XIA
		// input file name
		TString input_file_name;
		input_file_name.Form(
			"%s%stafd%d-map-%04u.root",
			kGenerateDataPath,
			kMappingDir,
			taf_index,
			alpha_calibration_run[taf_index]
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
		// input map event
		DssdMapEvent event;
		// setup input branches
		event.SetupInput(ipt);

		// total number of entries
		long long entries = ipt->GetEntries();
		// 1/100 of entries
		long long entry100 = entries / 100 + 1;
		// show start
		printf("Filling histogram   0%%");
		fflush(stdout);
		// fill energy to histogram
		for (long long entry = 0; entry < entries; ++entry) {
			// show process
			if (entry % entry100 == 0) {
				printf("\b\b\b\b%3lld%%", entry / entry100);
				fflush(stdout);
			}
			ipt->GetEntry(entry);
			// jump if back side
			if (event.side == 1) continue;
			// jump if energy out of range
			if (event.energy < 10000 || event.energy > 20000) continue;
			hist_channel[event.strip].Fill(event.energy);
			hist_energy[event.strip].Fill(
				cali_p0[event.strip] + cali_p1[event.strip] * event.energy
			);

			// trick to prevent error in bad strip
			if (taf_index == 5 && event.strip == 14) {
				hist_channel[15].Fill(event.energy);
			}
		}
		// show finish
		printf("\b\b\b\b100%%\n");
		// close input file
		ipf.Close();
	}

	TCanvas *c1 = new TCanvas("c1", "c1", 1920, 1080);
	c1->cd();
	TString pdf_name = TString::Format(
		"%simages/tafd%d-alpha-peaks.pdf",
		kGenerateDataPath,
		taf_index
	);
	c1->Print(pdf_name+"[");
	for (int i = 0; i < 16; ++i) {
		hist_channel[i].Draw();
		c1->Print(pdf_name);
	}
	for (int i = 0; i < 16; ++i) {
		hist_energy[i].Draw();
		c1->Print(pdf_name);
	}
	c1->Print(pdf_name+"]");

	opf.cd();
	for (int i = 0; i < 16; ++i) hist_channel[i].Write();
	for (int i = 0; i < 16; ++i) hist_energy[i].Write();
	opf.Close();
	return 0;
}