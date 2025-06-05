// spectrum 程序版本3，仅从 channel_v2 读取
#include <iostream>
#include <fstream>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TF1.h>
#include <Math/Vector3D.h>
#include <TCanvas.h>
#include <TPaveStats.h>

#include "include/event/channel_v2_event.h"
#include "include/ppac_track.h"

using namespace ribll;

constexpr double be_excited_energy[4] = {
	0.0, 3.368, 6.179, 7.542
};

constexpr double correct_tafd_thickness[6] = {
	166.0, 160.0, 150.0, 160.0, 164.0, 170.0
};

// Q value correct
constexpr double q_tafd_correct[6] = {
	0.3, -0.3, -0.4, 0.0, 0.4, 0.4
};
constexpr double q_csi_correct[12] = {
	-0.46, -0.25, 0.18, 0.49,
	0.49, 0.55, 0.28, 0.30,
	-0.04, -0.10, -0.12, -0.64
};

constexpr double exmin = 12.1;
constexpr double exmax = 27.1;

int main(int argc, char **argv) {
	if (argc > 1) {
		std::cout << "Usage: " << argv[0] << "\n";
		return -1;
	}

	// input file name
	TString input_file_name = TString::Format(
		"%s%sC14-10Be-4He-2H-v2.root",
		kGenerateDataPath,
		kChannelDir
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
	ChannelV2Event channel;
	// setup input branches
	channel.SetupInput(ipt);

	// output file name
	TString output_file_name = TString::Format(
		"%s%sq-slice-excitation.root",
		kGenerateDataPath,
		kSpectrumDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// excitation spectrum for differnece slices
	TH1F ex0[6] = {
		TH1F("x0q0", "ex0 for q from -10 to -11", 50, exmin, exmax),
		TH1F("x0q1", "ex0 for q from -11 to -12", 50, exmin, exmax),
		TH1F("x0q2", "ex0 for q from -12 to -13", 50, exmin, exmax),
		TH1F("x0q3", "ex0 for q from -13 to -14", 50, exmin, exmax),
		TH1F("x0q4", "ex0 for q from -14 to -15", 50, exmin, exmax),
		TH1F("x0q5", "ex0 for q from -15 to -16", 50, exmin, exmax),
	};
	TH1F ex1[8] = {
		TH1F("x1q0", "ex1 for q from -14 to -14.5", 50, exmin, exmax),
		TH1F("x1q1", "ex1 for q from -14.5 to -15", 50, exmin, exmax),
		TH1F("x1q2", "ex1 for q from -15 to -15.5", 50, exmin, exmax),
		TH1F("x1q3", "ex1 for q from -15.5 to -16", 50, exmin, exmax),
		TH1F("x1q4", "ex1 for q from -16 to -16.5", 50, exmin, exmax),
		TH1F("x1q5", "ex1 for q from -16.5 to -17", 50, exmin, exmax),
		TH1F("x1q6", "ex1 for q from -17 to -17.5", 50, exmin, exmax),
		TH1F("x1q7", "ex1 for q from -17.5 to -18", 50, exmin, exmax),
	};
	TH1F ex2[4] = {
		TH1F("x2q0", "ex2 for q from -16 to -17", 50, exmin, exmax),
		TH1F("x2q1", "ex2 for q from -17 to -18", 50, exmin, exmax),
		TH1F("x2q2", "ex2 for q from -18 to -19", 50, exmin, exmax),
		TH1F("x2q3", "ex2 for q from -19 to -20", 50, exmin, exmax),
	};
	// output tree
	TTree opt("tree", "spectrum v3");
	// output data
	int valid, ppac_flag, target_flag, taf_flag;
	// kinetic energy with target energy lost
	double be_kinetic, he_kinetic, d_kinetic;
	double c_momentum, c_kinetic;
	// Q value
	double q;
	// 10Be state, under different T0D2 energy method
	int be_state;
	// excited energy with target energy lost
	double excited_energy;

	// setup output branches
	opt.Branch("valid", &valid, "valid/I");
	opt.Branch("ppac_flag", &ppac_flag, "pflag/I");
	opt.Branch("taf_flag", &taf_flag, "tflag/I");
	opt.Branch("target_flag", &target_flag, "tarflag/I");
	opt.Branch("be_kinetic", &be_kinetic, "bek/D");
	opt.Branch("he_kinetic", &he_kinetic, "hek/D");
	opt.Branch("d_kinetic", &d_kinetic, "dk/D");
	opt.Branch("c_momentum", &c_momentum, "cp/D");
	opt.Branch("c_kinetic", &c_kinetic, "ck/D");
	opt.Branch("q", &q, "q/D");
	opt.Branch("be_state", &be_state, "bes/I");
	opt.Branch("excited_energy", &excited_energy, "ex/D");

	// loop to process
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		// get data
		ipt->GetEntry(entry);

		if (channel.hole != 0) continue;
		// get flags
		valid = 0;
		if (!channel.t0_valid) valid |= 1;
		taf_flag = 0;
		if (!channel.tafcsi_valid && channel.tafd_front_strip >= 13) {
			taf_flag = 1;
		} else if (channel.tafcsi_valid && channel.tafd_front_strip <= 13) {
			taf_flag = 2;
		} else {
			valid |= 2;
		}

		if (!channel.ppac_valid) valid |= 4;
		if (pow(channel.tx-2.5, 2.0) + pow(channel.ty+2.0, 2.0) < 196.0) {
			target_flag = 1;
		} else {
			target_flag = 0;
			valid |= 8;
		}

		double &tx = channel.tx;
		double &ty = channel.ty;

		// 10Be direction vector
		ROOT::Math::XYZVector d_be(
			channel.fragment_x[0] - tx,
			channel.fragment_y[0] - ty,
			100.0
		);
		d_be = d_be.Unit();
		// 4He direction vector
		ROOT::Math::XYZVector d_he(
			channel.fragment_x[1] - tx,
			channel.fragment_y[1] - ty,
			100.0
		);
		d_he = d_he.Unit();
		// 2H direction vector
		ROOT::Math::XYZVector d_d(
			channel.recoil_x - tx,
			channel.recoil_y - ty,
			135.0
		);
		d_d = d_d.Unit();

		// sum up
		be_kinetic = channel.fragment_kinetic[0];
		// 4He kinetic energy
		he_kinetic = channel.fragment_kinetic[1];
		// deutron kinetic energy
		d_kinetic = channel.recoil_kinetic;

		// 10Be momentum
		double be_momentum = MomentumFromKinetic(mass_10be, be_kinetic);
		ROOT::Math::XYZVector p_be = d_be * be_momentum;

		// 4He momentum
		double he_momentum = MomentumFromKinetic(mass_4he, he_kinetic);
		ROOT::Math::XYZVector p_he = d_he * he_momentum;

		// 2H momentum
		double d_momentum = MomentumFromKinetic(mass_2h, d_kinetic);
		ROOT::Math::XYZVector p_d = d_d * d_momentum;

		// beam 14C momentum vector
		ROOT::Math::XYZVector p_c = p_be + p_he + p_d;

		// 14C momentum
		c_momentum = p_c.R();
		// 14C kinematic energy
		c_kinetic =
			sqrt(pow(c_momentum, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

		// three body Q value
		q = be_kinetic + he_kinetic + d_kinetic - c_kinetic;

		if (c_kinetic < 360.0) valid |= 16;

		// correct Q
		if (taf_flag == 1) {
			q += q_tafd_correct[channel.taf_index];
		} else if (taf_flag == 2) {
			q -= q_csi_correct[channel.csi_index];
		}

		if (valid != 0) continue;

		// state 0
		if (q > -16 && q < -10) {
			be_momentum = MomentumFromKinetic(mass_10be, be_kinetic);
			p_be = d_be * be_momentum;
			// excited 14C momentum vector
			ROOT::Math::XYZVector p_excited_c = p_be + p_he;
			// excited 14C momentum
			double excited_c_momentum = p_excited_c.R();
			// excited 14C total energy
			double excited_c_energy = (be_kinetic + mass_10be)
				+ (he_kinetic + mass_4he);
			// excited 14C mass
			double excited_c_mass = sqrt(
				pow(excited_c_energy, 2.0) - pow(excited_c_momentum, 2.0)
			);
			// excited energy of 14C
			excited_energy = excited_c_mass - mass_14c;
			int slice = 5 - int(q + 16);
			ex0[slice].Fill(excited_energy);
			be_state = 0;
			opt.Fill();
		}

		// state 1
		if (q > -18 && q < -14) {
			be_momentum = MomentumFromKinetic(mass_10be, be_kinetic);
			p_be = d_be * be_momentum;
			// excited 14C momentum vector
			ROOT::Math::XYZVector p_excited_c = p_be + p_he;
			// excited 14C momentum
			double excited_c_momentum = p_excited_c.R();
			// excited 14C total energy
			double excited_c_energy = (be_kinetic + mass_10be + be_excited_energy[1])
				+ (he_kinetic + mass_4he);
			// excited 14C mass
			double excited_c_mass = sqrt(
				pow(excited_c_energy, 2.0) - pow(excited_c_momentum, 2.0)
			);
			// excited energy of 14C
			excited_energy = excited_c_mass - mass_14c;
			int slice = 7 - int((q + 18) * 2.0);
			ex1[slice].Fill(excited_energy);
			be_state = 1;
			opt.Fill();
		}


		// state 2
		if (q > -20 && q < -16) {
			be_momentum = MomentumFromKinetic(mass_10be, be_kinetic);
			p_be = d_be * be_momentum;
			// excited 14C momentum vector
			ROOT::Math::XYZVector p_excited_c = p_be + p_he;
			// excited 14C momentum
			double excited_c_momentum = p_excited_c.R();
			// excited 14C total energy
			double excited_c_energy = (be_kinetic + mass_10be + be_excited_energy[2])
				+ (he_kinetic + mass_4he);
			// excited 14C mass
			double excited_c_mass = sqrt(
				pow(excited_c_energy, 2.0) - pow(excited_c_momentum, 2.0)
			);
			// excited energy of 14C
			excited_energy = excited_c_mass - mass_14c;
			int slice = 3 - int(q + 20);
			ex2[slice].Fill(excited_energy);
			be_state = 2;
			opt.Fill();
		}
	}

	// save
	opf.cd();
	opt.Write();
	for (size_t i = 0; i < 6; ++i) ex0[i].Write();
	for (size_t i = 0; i < 8; ++i) ex1[i].Write();
	for (size_t i = 0; i < 4; ++i) ex2[i].Write();

	TCanvas *c0 = new TCanvas("c0");
	TPad *pads0[6];
	for (int i = 0; i < 6; ++i) {
		pads0[i] = new TPad(
			TString::Format("pad0%d", i), "",
			0, (16-i*3)/20.0, 1, (19-i*3)/20.0
		);
		pads0[i]->Draw();
	}
	for (int i = 0; i < 6; ++i) {
		pads0[i]->SetTopMargin(0);
		pads0[i]->SetBottomMargin(0);
	}
	pads0[5]->SetBottomMargin(0.01);
	for (int i = 0; i < 6; ++i) {
		pads0[i]->cd();
		ex0[i].SetTitle("");
		ex0[i].Draw();
		pads0[i]->Update();
		TPaveStats *stats = (TPaveStats*)ex0[i].FindObject("stats");
		if (stats) {
			stats->SetX1NDC(0.75);
			stats->SetX2NDC(0.89);
			stats->SetY1NDC(0.70);
			stats->SetY2NDC(0.98);
			stats->Draw();
		}
	}
	c0->Write("c0");

	TCanvas *c1 = new TCanvas("c1");
	TPad *pads1[8];
	for (int i = 0; i < 8; ++i) {
		pads1[i] = new TPad(
			TString::Format("pad1%d", i), "",
			0, (22-i*3)/26.0, 1, (25-i*3)/26.0
		);
		pads1[i]->Draw();
	}
	for (int i = 0; i < 8; ++i) {
		pads1[i]->SetTopMargin(0);
		pads1[i]->SetBottomMargin(0);
	}
	pads1[7]->SetBottomMargin(0.01);
	for (int i = 0; i < 8; ++i) {
		pads1[i]->cd();
		ex1[i].SetTitle("");
		ex1[i].Draw();
		pads1[i]->Update();
		TPaveStats *stats = (TPaveStats*)ex1[i].FindObject("stats");
		if (stats) {
			stats->SetX1NDC(0.75);
			stats->SetX2NDC(0.89);
			stats->SetY1NDC(0.70);
			stats->SetY2NDC(0.98);
			stats->Draw();
		}
	}
	c1->Write("c1");

	TCanvas *c2 = new TCanvas("c2");
	TPad *pads2[4];
	for (int i = 0; i < 4; ++i) {
		pads2[i] = new TPad(
			TString::Format("pad2%d", i), "",
			0, (10-i*3)/14.0, 1, (13-i*3)/14.0
		);
		pads2[i]->Draw();
	}
	for (int i = 0; i < 4; ++i) {
		pads2[i]->SetTopMargin(0);
		pads2[i]->SetBottomMargin(0);
	}
	pads2[3]->SetBottomMargin(0.01);
	for (int i = 0; i < 4; ++i) {
		pads2[i]->cd();
		ex2[i].SetTitle("");
		ex2[i].Draw();
		pads2[i]->Update();
		TPaveStats *stats = (TPaveStats*)ex2[i].FindObject("stats");
		if (stats) {
			stats->SetX1NDC(0.75);
			stats->SetX2NDC(0.89);
			stats->SetY1NDC(0.70);
			stats->SetY2NDC(0.98);
			stats->Draw();
		}
	}
	c2->Write("c2");

	// close files
	opf.Close();
	ipf.Close();

	return 0;
}