#include <iostream>

#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TString.h>
#include <Math/Vector3D.h>

#include "include/event/generate_event.h"
#include "include/event/dssd_event.h"
#include "include/event/csi_event.h"
#include "include/calculator/range_energy_calculator.h"

using namespace ribll;

const std::pair<double, double> tafd_phi_ranges[6] = {
	{117.6*TMath::DegToRad(), 62.4*TMath::DegToRad()},
	{57.6*TMath::DegToRad(), 2.4*TMath::DegToRad()},
	{-2.4*TMath::DegToRad(), -57.6*TMath::DegToRad()},
	{-62.4*TMath::DegToRad(), -117.6*TMath::DegToRad()},
	{-122.4*TMath::DegToRad(), -177.6*TMath::DegToRad()},
	{177.6*TMath::DegToRad(), 122.4*TMath::DegToRad()}
};
constexpr double tafd_threshold[6] = {0.6, 0.5, 0.6, 0.6, 0.5, 0.5};
constexpr double csi_threshold[12] = {
	1200.0, 1200.0,
	1200.0, 1200.0,
	1200.0, 1400.0,
	1000.0, 1200.0,
	1100.0, 1200.0,
	1100.0, 1200.0
};

int main(int argc, char **argv) {
	int run = 0;
	if (argc > 1) {
		run = atoi(argv[1]);
	}
	if (run < 0 || run > 2) {
		std::cout << "Usage: " << argv[0] << "[run]\n"
			<< "  run        run number, default is 0\n";
	}

	// input generate data file name
	TString generate_file_name = TString::Format(
		"%s%sgenerate-%04d.root",
		kGenerateDataPath,
		kSimulateDir,
		run
	);
	// generate data file
	TFile generate_file(generate_file_name, "read");
	// input generate tree
	TTree *generate_tree = (TTree*)generate_file.Get("tree");
	if (!generate_tree) {
		std::cerr << "Error: Get tree from "
			<< generate_file_name << " failed.\n";
		return -1;
	}
	// generate event
	GenerateEvent event;
	// setup input branches
	event.SetupInput(generate_tree);

	// output TAFD files name
	TString tafd_file_names[6];
	for (int i = 0; i < 6; ++i) {
	 	tafd_file_names[i].Form(
			"%s%stafd%d-merge-sim-ta-%04d.root",
			kGenerateDataPath,
			kMergeDir,
			i,
			run
		);
	}
	// output TAFD files
	TFile *tafd_files[6];
	for (int i = 0; i < 6; ++i) {
		tafd_files[i] = new TFile(tafd_file_names[i], "recreate");
	}
	// output TAFD trees
	TTree *tafd_trees[6];
	for (int i = 0; i < 6; ++i) {
		tafd_files[i]->cd();
		tafd_trees[i] = new TTree("tree", "simulated TAFD data");
	}
	// output TAFD merge events
	AdssdMergeEvent tafd_events[6];
	// setup output branches
	for (int i = 0; i < 6; ++i) {
		tafd_events[i].SetupOutput(tafd_trees[i]);
	}

	// output TAFCsI file name
	TString csi_file_name = TString::Format(
		"%s%stafcsi-fundamental-sim-ta-%04d.root",
		kGenerateDataPath,
		kFundamentalDir,
		run
	);
	// output TAFCsI file
	TFile csi_file(csi_file_name, "recreate");
	// output TAFCsI tree
	TTree csi_tree("tree", "simulated TAFCsI data");
	// output TAFCsI event
	CircularCsiFundamentalEvent csi_event;
	// setup output branches
	csi_event.SetupOutput(&csi_tree);

	// initialize random number generator
	TRandom3 generator(0);

	// initialize calculators
	elc::RangeEnergyCalculator h2_calculator("2H", "Si");


	// initialize useless variables
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 4; ++j) {
			tafd_events[i].time[j] = 0.0;
			tafd_events[i].decode_entry[j] = 0;
		}
	}
	csi_event.cfd_flag = 0;
	for (int i = 0; i < 12; ++i) {
		csi_event.decode_entry[i] = 0;
	}

	// show start
	printf("Simulating TAF detect   0%%");
	// total entries, fow showing process
	long long entries = generate_tree->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100;
	// detecting simulation
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get entry
		generate_tree->GetEntry(entry);

		// initialize
		for (int i = 0; i < 6; ++i) {
			tafd_events[i].hit = 0;
		}
		csi_event.match = false;
		for (int i = 0; i < 12; ++i) {
			csi_event.time[i] = -1e5;
			csi_event.energy[i] = 0.0;
		}

		// TAF position
		ROOT::Math::XYZVector taf_position(
			event.rx, event.ry, event.rz
		);
		// TAF index
		int taf_index = int(8 - event.recoil_phi / pi * 3.0) % 6;
		// middle phi
		double mid_phi = (
			tafd_phi_ranges[taf_index].second
			+ tafd_phi_ranges[taf_index].first
		) / 2.0;
		// circle center of ADSSD
		ROOT::Math::XYZVector center(
			34.4*cos(mid_phi), 34.4*sin(mid_phi), 135.0
		);
		// position relate to cirlce center
		auto relative_position = taf_position - center;
		// radius from cirlce center
		double radius = relative_position.R();
		// phi value
		double phi = relative_position.Phi();
		// out of range
		if (
			radius < 32.6
			|| radius > 135.1
			|| phi < tafd_phi_ranges[taf_index].second
			|| phi > tafd_phi_ranges[taf_index].first
		) {
			for (int i = 0; i < 6; ++i) tafd_trees[i]->Fill();
			csi_tree.Fill();
			continue;
		}
		// cover range of single ring strip
		double tafd_ring_width = (135.1 - 32.6) / 16.0;
		// TAFD ring strip
		int tafd_ring_strip = int((radius-32.6) / tafd_ring_width);
		// cover range of single phi strip
		double tafd_phi_width = -55.2 * TMath::DegToRad() / 8.0;
		// TAFD phi strip
		int tafd_phi_strip = int(
			(phi - tafd_phi_ranges[taf_index].first) / tafd_phi_width
		);

		// convert strips to position back
		double merge_radius =
			tafd_ring_width * (tafd_ring_strip + 0.5) + 32.6;
		double merge_phi =
			tafd_phi_width * (tafd_phi_strip + 0.5)
			+ tafd_phi_ranges[taf_index].first;
		ROOT::Math::XYZVector merge_position(
			merge_radius * cos(merge_phi),
			merge_radius * sin(merge_phi),
			0.0
		);
		merge_position += center;

		// TAFCsI index
		int csi_index = 2 * taf_index;
		if (tafd_phi_strip >= 4) csi_index += 1;
		// calculate taf energy
		double recoil_range =
			h2_calculator.Range(event.recoil_kinetic_after_target);
		double thickness = tafd_thickness[taf_index];
		int taf_layer;
		double taf_lost_energy[2];
		if (recoil_range < thickness / cos(event.recoil_theta)) {
			taf_layer = 0;
			taf_lost_energy[0] = event.recoil_kinetic_after_target;
			taf_lost_energy[1] = 0.0;
		} else {
			taf_layer = 1;
			taf_lost_energy[1] = h2_calculator.Energy(
				recoil_range - thickness / cos(event.recoil_theta)
			);
			taf_lost_energy[0] =
				event.recoil_kinetic_after_target - taf_lost_energy[1];
		}

		// consider energy resolution
		double tafd_energy = taf_lost_energy[0] + generator.Gaus(0.0, 0.05);
		double csi_energy = taf_lost_energy[1] + generator.Gaus(0.0, 1.0);
		// double tafd_energy = taf_lost_energy[0];
		// double csi_energy = taf_lost_energy[1];

		// convert energy to channel
		double csi_channel =
			power_csi_param[csi_index][0]
			* pow(csi_energy, power_csi_param[csi_index][1])
			+ power_csi_param[csi_index][2];

		// fill event
		// fill TAFD merge event
		if (tafd_energy > tafd_threshold[taf_index]) {
			tafd_events[taf_index].hit = 1;
			tafd_events[taf_index].radius[0] = merge_position.R();
			tafd_events[taf_index].theta[0] = merge_position.Theta();
			tafd_events[taf_index].phi[0] = merge_position.Phi();
			tafd_events[taf_index].front_strip[0] = tafd_ring_strip;
			tafd_events[taf_index].back_strip[0] = tafd_phi_strip;
			tafd_events[taf_index].energy[0] = tafd_energy;
// std::cout << "Entry " << entry << ", layer " << taf_layer << ": Generate energy " << event.recoil_kinetic_after_target <<  ", range "
// 	<< recoil_range << ", taf index " << taf_index << ", csi_index " << csi_index << ", thick " << thickness
// 	<< ", effect thickness " << thickness / cos(event.recoil_theta) << ", tafd energy " << tafd_energy << ", csi energy " << csi_energy
// 	<< ", csi channel " << csi_channel << "\n";

		}
		// fill TAFCsI event
		if (taf_layer == 1 && csi_channel > csi_threshold[csi_index]) {
			csi_event.match = true;
			csi_event.time[csi_index] = 0.0;
			csi_event.energy[csi_index] = csi_channel;
// std::cout << "Entry " << entry << ", layer " << taf_layer << ": Generate energy " << event.recoil_kinetic_after_target <<  ", range "
// 	<< recoil_range << ", taf index " << taf_index << ", csi_index " << csi_index << ", thick " << thickness
// 	<< ", effect thickness " << thickness / cos(event.recoil_theta)  << ", tafd energy " << tafd_energy << ", csi energy " << csi_energy
// 	<< ", csi channel " << csi_channel << "\n";
		}

		// fill to tree
		for (int i = 0; i < 6; ++i) tafd_trees[i]->Fill();
		csi_tree.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");


	// save trees and close output files
	for (int i = 0; i < 6; ++i) {
		tafd_files[i]->cd();
		tafd_trees[i]->Write();
		tafd_files[i]->Close();
	}
	csi_file.cd();
	csi_tree.Write();
	csi_file.Close();
	// close input file
	generate_file.Close();

	return 0;
}