#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>

#include "include/event/generate_event.h"
#include "include/event/threebody_info_event.h"

using namespace ribll;

int main() {
	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody-2H-sim-2-2.root",
		kGenerateDataPath, kSpectrumDir
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
	int be_state[4];
	double be_kinetic_target[3], he_kinetic_target[3];
	double excited_energy_target[3];
	// setup input branches
	ipt->SetBranchAddress("be_state", be_state);
	ipt->SetBranchAddress("be_kinetic_target", be_kinetic_target);
	ipt->SetBranchAddress("he_kinetic_target", he_kinetic_target);
	ipt->SetBranchAddress("excited_energy_target", excited_energy_target);

	// generate file
	TString generate_file_name = TString::Format(
		"%s%sgenerate-0002.root", kGenerateDataPath, kSimulateDir
	);
	// add friend
	ipt->AddFriend("g=tree", generate_file_name);
	// input data
	double generate_kinetic[2];
	double c14ex;
	// setup input branches
	ipt->SetBranchAddress("g.fragment_kinetic_in_target", generate_kinetic);
	ipt->SetBranchAddress("g.c14_excited", &c14ex);

	// info file
	TString info_file_name = TString::Format(
		"%s%sthreebody-2H-sim-2.root", kGenerateDataPath, kInformationDir
	);
	// add friend
	ipt->AddFriend("info=tree", info_file_name);
	// input event
	ThreeBodyInfoEvent info;
	// setup input branches
	info.SetupInput(ipt, "info.");

	// output file name
	TString output_file_name = TString::Format(
		"%s%sd2_calculate_resolution.root", kGenerateDataPath, kShowDir
	);
	// ouptut file
	TFile opf(output_file_name, "recreate");

	// histogram of difference between generate and meansured
	// 10Be stop in T0D2
	TH1F h_d2_bek0("hbek0l1", "10Be stop in D2, fk[0]-bek[0]", 100, -10, 10);
	// 10Be stop in T0D3
	TH1F h_d3_bek0("hbek0l2", "10Be stop in D3, fk[0]-bek[0]", 100, -10, 10);
	// 4He stop in T0D2
	TH1F h_d2_hek0("hhek0l1", "4He stop in D2, fk[1]-hek[0]", 100, -10, 10);
	// 4He stop in T0D3
	TH1F h_d3_hek0("hhek0l2", "4He stop in D3, fk[1]-hek[0]", 100, -10, 10);
	// 4He stop in SSD
	TH1F h_ssd_hek0("hhek0l3", "4He stop in SSD, fk[1]-hek[0]", 100, -10, 10);
	// hitogram pointers
	TH1F *hist_k0[5] = {
		&h_d2_bek0, &h_d3_bek0, &h_d2_hek0, &h_d3_hek0, &h_ssd_hek0
	};

	// histogram of difference between generate and calculated
	// 10Be stop in T0D2
	TH1F h_d2_bek2("hbek2l1", "10Be stop in D2, fk[0]-bek[2]", 100, -10, 10);
	// 10Be stop in T0D3
	TH1F h_d3_bek2("hbek2l2", "10Be stop in D3, fk[0]-bek[2]", 100, -10, 10);
	// 4He stop in T0D2
	TH1F h_d2_hek2("hhek2l1", "4He stop in D2, fk[1]-hek[2]", 100, -10, 10);
	// 4He stop in T0D3
	TH1F h_d3_hek2("hhek2l2", "4He stop in D3, fk[1]-hek[2]", 100, -10, 10);
	// 4He stop in SSD
	TH1F h_ssd_hek2("hhek2l3", "4He stop in SSD, fk[1]-hek[2]", 100, -10, 10);
	// histogram pointers
	TH1F *hist_k2[5] = {
		&h_d2_bek2, &h_d3_bek2, &h_d2_hek2, &h_d3_hek2, &h_ssd_hek2
	};

	// histogram of difference between measured and calculated
	// 10Be stop in T0D2
	TH1F h_d2_bekd("hbekdl1", "10Be stop in D2, bek[0]-bek[2]", 100, -10, 10);
	// 10Be stop in T0D3
	TH1F h_d3_bekd("hbekdl2", "10Be stop in D3, bek[0]-bek[2]", 100, -10, 10);
	// 4He stop in T0D2
	TH1F h_d2_hekd("hhekdl1", "4He stop in D2, hek[0]-hek[2]", 100, -10, 10);
	// 4He stop in T0D3
	TH1F h_d3_hekd("hhekdl2", "4He stop in D3, hek[0]-hek[2]", 100, -10, 10);
	// 4He stop in SSD
	TH1F h_ssd_hekd("hhekdl3", "4He stop in SSD, hek[0]-hek[2]", 100, -10, 10);
	// histogram pointers
	TH1F *hist_kd[5] = {
		&h_d2_bekd, &h_d3_bekd, &h_d2_hekd, &h_d3_hekd, &h_ssd_hekd
	};

	// histogram of difference between generated and measured
	TH1F h_ex0("hex0", "c14ex-ex[0]", 100, -10, 10);
	TH1F h_ex2("hex2", "c14ex-ex[2]", 100, -10, 10);


	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling hisotgrams   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get data
		ipt->GetEntry(entry);

		if ((info.ppac_flag & 1) != 1) continue;
		if (info.taf_flag != 0) continue;
		if (info.layer[0] == 1) {
			h_d2_bek0.Fill(generate_kinetic[0] - be_kinetic_target[0]);
			h_d2_bek2.Fill(generate_kinetic[0] - be_kinetic_target[2]);
			h_d2_bekd.Fill(be_kinetic_target[0] - be_kinetic_target[2]);
		} else if (info.layer[0] == 2) {
			h_d3_bek0.Fill(generate_kinetic[0] - be_kinetic_target[0]);
			h_d3_bek2.Fill(generate_kinetic[0] - be_kinetic_target[2]);
			h_d3_bekd.Fill(be_kinetic_target[0] - be_kinetic_target[2]);
		}
		if (info.layer[1] == 1) {
			h_d2_hek0.Fill(generate_kinetic[1] - he_kinetic_target[0]);
			h_d2_hek2.Fill(generate_kinetic[1] - he_kinetic_target[2]);
			h_d2_hekd.Fill(he_kinetic_target[0] - he_kinetic_target[2]);
		} else if (info.layer[1] == 2) {
			h_d3_hek0.Fill(generate_kinetic[1] - he_kinetic_target[0]);
			h_d3_hek2.Fill(generate_kinetic[1] - he_kinetic_target[2]);
			h_d3_hekd.Fill(he_kinetic_target[0] - he_kinetic_target[2]);
		} else if (info.layer[1] > 2) {
			h_ssd_hek0.Fill(generate_kinetic[1] - he_kinetic_target[0]);
			h_ssd_hek2.Fill(generate_kinetic[1] - he_kinetic_target[2]);
			h_ssd_hekd.Fill(he_kinetic_target[0] - he_kinetic_target[2]);
		}

		h_ex0.Fill(c14ex - excited_energy_target[0]);
		h_ex2.Fill(c14ex - excited_energy_target[0]);
	}
	printf("\b\b\b\b100%%\n");

	double sigma_k0[5], mean_k0[5];
	// fit
	for (int i = 0; i < 5; ++i) {
		TF1 *f1 = new TF1(TString::Format("f0_%d", i), "gaus", -5, 5);
		hist_k0[i]->Fit(f1, "RQ+");
		mean_k0[i] = f1->GetParameter(1);
		sigma_k0[i] = f1->GetParameter(2);
	}
	// show fitting result
	std::string title_k0[5] = {
		"D2-fk[0]-bek[0]", "D3-fk[0]-bek[0]",
		"D2-fk[1]-hek[0]", "D3-fk[1]-hek[0]", "SSD-fk[1]-hek[0]"
	};
	for (int i = 0; i < 5; ++i) {
		std::cout << title_k0[i] << ", "
			<< mean_k0[i] << ", "
			<< sigma_k0[i] << "\n";
	}

	double sigma_k2[5], mean_k2[5];
	// fit
	for (int i = 0; i < 5; ++i) {
		TF1 *f1 = new TF1(TString::Format("f2_%d", i), "gaus", -5, 5);
		hist_k2[i]->Fit(f1, "RQ+");
		mean_k2[i] = f1->GetParameter(1);
		sigma_k2[i] = f1->GetParameter(2);
	}
	// show fitting result
	std::string title_k2[5] = {
		"D2-fk[0]-bek[2]", "D3-fk[0]-bek[2]",
		"D2-fk[1]-hek[2]", "D3-fk[1]-hek[2]", "SSD-fk[1]-hek[2]"
	};
	for (int i = 0; i < 5; ++i) {
		std::cout << title_k2[i] << ", "
			<< mean_k2[i] << ", "
			<< sigma_k2[i] << "\n";
	}

	double sigma_kd[5], mean_kd[5];
	// fit
	for (int i = 0; i < 5; ++i) {
		TF1 *f1 = new TF1(TString::Format("fd_%d", i), "gaus", -5, 5);
		hist_kd[i]->Fit(f1, "RQ+");
		mean_kd[i] = f1->GetParameter(1);
		sigma_kd[i] = f1->GetParameter(2);
	}
	// show fitting result
	std::string title_kd[5] = {
		"D2-bek[0]-bek[2]", "D3-bek[0]-bek[2]",
		"D2-hek[0]-hek[2]", "D3-hek[0]-hek[2]", "SSD-hek[0]-hek[2]"
	};
	for (int i = 0; i < 5; ++i) {
		std::cout << title_kd[i] << ", "
			<< mean_kd[i] << ", "
			<< sigma_kd[i] << "\n";
	}



	for (int i = 0; i < 5; ++i) hist_k0[i]->Write();
	for (int i = 0; i < 5; ++i) hist_k2[i]->Write();
	for (int i = 0; i < 5; ++i) hist_kd[i]->Write();
	h_ex0.Write();
	h_ex2.Write();
	// close file
	opf.Close();
	ipf.Close();
	return 0;
}