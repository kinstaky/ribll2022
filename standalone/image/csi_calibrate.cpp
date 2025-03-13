#include <iostream>

#include <TGraph.h>
#include <TFile.h>

#include "include/calculator/csi_energy_calculator.h"
#include "include/event/ta_event.h"

using namespace ribll;

void PrintUsage(const char* name) {
	std::cout << "Usage: " << name << " taf_index\n";
}

int main(int argc, char **argv) {
	if (argc != 2) {
		PrintUsage(argv[0]);
        return -1;
	}
	int taf_index = atoi(argv[1]);

	// input file name
	TString input_file_name = TString::Format(
		"%s%staf%d-telescope-sim-ta-0002.root",
		kGenerateDataPath,
		kTelescopeDir,
		taf_index
	);
	TFile ipf(input_file_name, "read");
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
            << input_file_name << " failed.\n";
        return -1;
	}
	TaEvent taf_event;
	taf_event.SetupInput(ipt);

	TString output_file_name = TString::Format(
		"%s%staf%dcsi-calibrate-sim-0002.root",
		kGenerateDataPath,
		kImageDir,
		taf_index
	);
	TFile opf(output_file_name, "recreate");
	TGraph gda, gdb;

	// calculator
	elc::CsiEnergyCalculator h2_csi_calculator("2H");

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling CsI calibration graph   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get events
		ipt->GetEntry(entry);
		if (taf_event.num != 1) continue;
		if (taf_event.flag[0] != 0x3 && taf_event.flag[0] != 0x5) continue;

		// CsI channel number
		double csi_channel = taf_event.energy[0][1];
		// CsI energy
		double csi_energy = h2_csi_calculator.Energy(
			taf_event.theta[0], taf_event.energy[0][0], tafd_thickness[taf_index]
		);
		if (taf_event.flag[0] == 0x3) {
			gda.AddPoint(csi_channel, csi_energy);
		} else {
			gdb.AddPoint(csi_channel, csi_energy);
		}
	}
	printf("\b\b\b\b100%%\n");

	opf.cd();
	gda.Write("gda");
	gdb.Write("gdb");

	opf.Close();
	ipf.Close();
}