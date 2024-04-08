#include <iostream>
#include <fstream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TF1.h>

#include "include/event/ta_event.h"
#include "include/event/particle_type_event.h"
#include "include/calculator/csi_energy_calculator.h"

using namespace ribll;


int main(int argc, char **argv) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " start_run end_run taf_index\n";
        return -1;
    }
    int start_run = atoi(argv[1]);
    int end_run = atoi(argv[2]);

	// output root file name
	TString output_file_name;
	output_file_name.Form(
		"%s%stafcsi-energy-range-%04d-%04d.root",
		kGenerateDataPath,
		kShowDir,
		start_run,
		end_run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output calibration graphs
	TGraph gdcali[12][40];
	std::vector<TH1F> hd_energy_range[12][40];
	for (int csi_index = 0; csi_index < 12; ++csi_index) {
		for (int t = 0; t < 40; ++t) {
			for (int i = 0; i < 48; ++i) {
				hd_energy_range[csi_index][t].emplace_back(
					TString::Format("he%dt%dr%d", csi_index, 131+t, i),
					TString::Format(
						"TAFD thickness %d, CsI energy in range %d-%d",
						131+t, 1400+200*i, 1600+200*i
					),
					600, 0, 60
				);
			}
		}
	}

	// CsI energy calculator
	elc::CsiEnergyCalculator h2_csi_calculator("2H");


	for (int taf_index = 0; taf_index < 6; ++taf_index) {
		// taf telescope event chain
		TChain taf_chain("taf", "chain of taf events");
		// particle type event chain
		TChain type_chain("type", "chain of particle type events");
		for (int run = start_run; run <= end_run; ++run) {
			if (run < 618) continue;
			if (run == 628) continue;
			if (run > 652 && run < 675) continue;
			if (run > 716 && run < 739) continue;
			if (run > 746) continue;

			taf_chain.AddFile(TString::Format(
				"%s%staf%d-telescope-ta-%04u.root/tree",
				kGenerateDataPath,
				kTelescopeDir,
				taf_index,
				run
			));
			type_chain.AddFile(TString::Format(
				"%s%staf%d-particle-type-ta-%04u.root/tree",
				kGenerateDataPath,
				kParticleIdentifyDir,
				taf_index,
				run
			));
		}
		taf_chain.AddFriend(&type_chain);
		// input TA telescope event
		TaEvent ta_event;
		// input particle type event
		ParticleTypeEvent type_event;
		// setup output branches
		ta_event.SetupInput(&taf_chain);
		type_event.SetupInput(&taf_chain, "type.");

		
		// total number of entries
		long long entries = taf_chain.GetEntries();
		// 1/100 of entries
		long long entry100 = entries / 100 + 1;
		// show start
		printf("Filling TAF%d CsI calibration graph   0%%", taf_index);
		fflush(stdout);
		for (long long entry = 0; entry < entries; ++entry) {
			// show process
			if (entry % entry100 == 0) {
				printf("\b\b\b\b%3lld%%", entry / entry100);
				fflush(stdout);
			}
			// get events
			taf_chain.GetEntry(entry);
			if (ta_event.num != 1) continue;
			if (ta_event.flag[0] != 0x3 && ta_event.flag[0] != 0x5) continue;
			// select only 2H
			if (type_event.mass[0] != 2 || type_event.charge[0] != 1) continue;

			// calculate the first index in graph array case by phi angle
			int csi_index = ta_event.flag[0] == 0x3 ? 0 : 1;
			csi_index += taf_index * 2;
			// CsI channel number
			double csi_channel = ta_event.energy[0][1];
			// loop thickness, from 131um to 170um, step 1um
			for (int thick = 0; thick < 40; ++thick) {
				double thickness = 131 + thick;
				// CsI energy
				double csi_energy = h2_csi_calculator.Energy(
					ta_event.theta[0], ta_event.energy[0][0], thickness
				);
				gdcali[csi_index][thick].AddPoint(csi_channel, csi_energy);
				if (csi_channel > 1400.0 && csi_channel < 11000.0) {
					int range_index = int((csi_channel - 1400.0) / 200.0);
					hd_energy_range[csi_index][thick][range_index].Fill(csi_energy);
				}
			}

		}
		// show finish
		printf("\b\b\b\b100%%\n");
	}

	double cali_param[40][12][2];

	// fit and calibrate
	for (int i = 0; i < 12; ++i) {
		for (int t = 0; t < 40; ++t) {
			TF1 *f1 = new TF1(
				TString::Format("f%dt%d", i, t+131), "pol1", 0.0, 11000.0
			);
			gdcali[i][t].Fit(f1, "RQ+");
			f1->GetParameters(cali_param[t][i]);
		}
	}

	// save calibration parameters
	// file name
	TString parameter_file_name = TString::Format(
		"%s%stafcsi-pol1-cali-param.txt", kGenerateDataPath, kCalibrationDir
	);
	// output stream
	std::ofstream fout(parameter_file_name.Data());
	if (!fout.good()) {
		std::cerr << "Error: Open parameter file "
			<< parameter_file_name << " failed.\n";
		return -1;
	}
	// write parameters
	for (int t = 0; t < 40; ++t) {
		for (int i = 0; i < 11; ++i) {
			fout << cali_param[t][i][0] << " "
				<< cali_param[t][i][1] << " ";
		}
		fout << cali_param[t][11][0] << " "
			<< cali_param[t][11][1] << "\n";
	}
	// close file
	fout.close();


	opf.cd();
	// save graphs
	for (int i = 0; i < 12; ++i) {
		for (int t = 0; t < 40; ++t) {
	        gdcali[i][t].Write(TString::Format("gd%dt%d", i, t+131));
		}
	}
	for (int i = 0; i < 12; ++i) {
		for (int t = 0; t < 40; ++t) {
			for (int j = 0; j < 48; ++j) {
				hd_energy_range[i][t][j].Write();
			}
		}
	}

	// close files
	opf.Close();

    return 0;
}