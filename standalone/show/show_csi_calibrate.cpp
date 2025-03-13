/*
 * 这个程序是用来看模拟或者实验数据中，通过 ADSSD 推导的 CsI 能量和 CsI 的通道的关系
 * 用来和 CsI 刻度的实验数据进行比较，以得到比较合理的 CsI 的能量分辨
 */
#include <iostream>

#include <TGraph.h>
#include <TFile.h>
#include <TChain.h>

#include "include/calculator/csi_energy_calculator.h"
#include "include/event/ta_event.h"
#include "include/event/particle_event.h"

using namespace ribll;


void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " [options] taf_index\n"
		"Options:\n"
		"  -h                Print this help information.\n"
		"  -s                Use simulated data.\n";
}


/// @brief parse arguments
/// @param[in] argc number of arguments
/// @param[in] argv arguments
/// @param[out] help need help
/// @param[out] sim use simulated data
/// @returns start index of positional arguments if succes, if failed returns
///		-argc (negative argc) for miss argument behind option,
/// 	or -index (negative index) for invalid arguemnt
///
int ParseArguments(
	int argc,
	char **argv,
	bool &help,
	bool &sim
) {
	// initialize
	help = false;
	sim = false;
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
		} else if (argv[result][1] == 's') {
			sim = true;
		} else {
			return -result;
		}
	}
	return result;
}


int ShowSimulate(int taf_index) {
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
		kShowDir,
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
			gda.AddPoint(csi_energy, csi_channel);
		} else {
			gdb.AddPoint(csi_energy, csi_channel);
		}
	}
	printf("\b\b\b\b100%%\n");

	opf.cd();
	gda.Write("gdpa");
	gdb.Write("gdpb");

	opf.Close();
	ipf.Close();

	return 0;
}


int Show(int taf_index) {
	TChain taf_chain("taf", "chain of taf events");
	TChain particle_chain("particle", "chain of particle events");
	for (unsigned int run = 618; run <= 746; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		if (run > 716 && run < 739) continue;

		taf_chain.AddFile(TString::Format(
			"%s%staf%d-telescope-ta-%04u.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			taf_index,
			run
		));
		particle_chain.AddFile(TString::Format(
			"%s%staf%d-particle-ta-v2-%04u.root/tree",
			kGenerateDataPath,
			kParticleDir,
			taf_index,
			run
		));
	}
	taf_chain.AddFriend(&particle_chain);
	// input taf event
	TaEvent taf_event;
	// input particle event
	ParticleEvent particle;
	// setup output branches
	taf_event.SetupInput(&taf_chain);
	particle.SetupInput(&taf_chain, "particle.");


	TString output_file_name = TString::Format(
		"%s%staf%dcsi-calibrate.root",
		kGenerateDataPath,
		kShowDir,
		taf_index
	);
	TFile opf(output_file_name, "recreate");
	TGraph gda, gdb;

	// calculator
	elc::CsiEnergyCalculator h2_csi_calculator("2H");

	// total number of entries
	long long entries = taf_chain.GetEntries();
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
		taf_chain.GetEntry(entry);
		if (particle.num != 1 || taf_event.num != 1) continue;
		if (taf_event.flag[0] != 0x3 && taf_event.flag[0] != 0x5) continue;
		if (particle.mass[0] != 2 || particle.charge[0] != 1) continue;

		// CsI channel number
		double csi_channel = taf_event.energy[0][1];
		// CsI energy
		double csi_energy = h2_csi_calculator.Energy(
			taf_event.theta[0], taf_event.energy[0][0], tafd_thickness[taf_index]
		);
		if (taf_event.flag[0] == 0x3) {
			gda.AddPoint(csi_energy, csi_channel);
		} else {
			gdb.AddPoint(csi_energy, csi_channel);
		}
	}
	printf("\b\b\b\b100%%\n");

	opf.cd();
	gda.Write("gdpa");
	gdb.Write("gdpb");

	opf.Close();

	return 0;
}


int main(int argc, char **argv) {
	if (argc > 3) {
		PrintUsage(argv[0]);
		return -1;
	}
	// help flag
	bool help = false;
	// simulated data
	bool sim = false;
	// parse arguments and get start index of positional arguments
	int pos_start = ParseArguments(argc, argv, help, sim);
	// need help
	if (help) {
		PrintUsage(argv[0]);
		return 0;
	}
	if (pos_start < 0) {
		if (-pos_start < argc) {
			std::cerr << "Error: Invaild option " << argv[-pos_start] << ".\n";
		} else {
			std::cerr << "Error: Option need parameter.\n";
		}
		PrintUsage(argv[0]);
		return -1;
	}
	if (pos_start >= argc) {
		// positional arguments less than 1
		std::cerr << "Error: Parameter case not found.\n";
		PrintUsage(argv[0]);
		return -1;
	}
	// TAF index
	int taf_index = atoi(argv[pos_start]);

	if (sim) {
		return ShowSimulate(taf_index);
	} else {
		return Show(taf_index);
	}

}