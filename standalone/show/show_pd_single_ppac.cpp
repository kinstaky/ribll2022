#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TH1F.h>
#include <Math/Vector3D.h>

#include "include/event/pd_info_event.h"
#include "include/ppac_track.h"

using namespace ribll;

double PdProcess(
	const PdInfoEvent &event,
	const double target_x,
	const double target_y,
	double &c14_kinetic,
	double &d_kinetic,
	double &c15_kinetic
) {
	// 14C kinetic energy
	// T0D1 energy
	c14_kinetic = t0_param[0][0] + t0_param[0][1] * event.c_channel[0];
	// T0D2 energy
	c14_kinetic += t0_param[1][0] + t0_param[1][1] * event.c_channel[1];

	double a0 = csi_param[event.csi_index][0];
	double a1 = csi_param[event.csi_index][1];
	double a2 = csi_param[event.csi_index][2];
	// 2H kinetic energy
	d_kinetic = event.tafd_energy + pow(
		(event.csi_channel - a2) / a0,
		1.0 / a1
	);

	// C14 momentum
	double c14_momentum = MomentumFromKinetic(mass_14c, c14_kinetic);
	// C14 momentum vector
	ROOT::Math::XYZVector p_c14(
		event.c_x[0] - target_x,
		event.c_y[0] - target_y,
		100.0
	);
	p_c14 = p_c14.Unit() * c14_momentum;

	// 2H momentum
	double d_momentum = MomentumFromKinetic(mass_2h, d_kinetic);
	// 2H momentum vector
	ROOT::Math::XYZVector p_d(
		event.d_x - target_x,
		event.d_y - target_y,
		135.0
	);
	p_d = p_d.Unit() * d_momentum;

	// beam 14C momentum vector
	ROOT::Math::XYZVector p_c15 = p_c14 + p_d;

	// 15C momentum
	double c15_momentum = p_c15.R();
	// 15C kinetic energy
	c15_kinetic =
		sqrt(pow(c15_momentum, 2.0) + pow(mass_15c, 2.0)) - mass_15c;

	// Q value
	double q = c14_kinetic + d_kinetic - c15_kinetic;

	return q;
}


int main() {
	// input file name
	TString input_file_name = TString::Format(
		"%s%spd.root", kGenerateDataPath, kInformationDir
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
	PdInfoEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// output file name
	TString output_file_name = TString::Format(
		"%s%spd-single-ppac.root", kGenerateDataPath, kShowDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// histograms of dX and dY of single and multiple tracked point
	TH1F *hist_dx[3], *hist_dy[3];
	for (int i = 0; i < 3; ++i) {
		hist_dx[i] = new TH1F(
			TString::Format("hdx%d", i), "dX", 100, -10, 10
		);
		hist_dy[i] = new TH1F(
			TString::Format("hdy%d", i), "dY", 100, -10, 10
		);
	}
	// histograms of Q value with multiple PPAC tracking
	TH1F hist_q_mpt(
		"hqmpt", "Q with multiple PPAC track", 150, -10, 5
	);
	// histograms of Q value with single PPAC tracking approximation (spta)
	TH1F *hist_q_spta[3];
	for (int i = 0; i < 3; ++i) {
		hist_q_spta[i] = new TH1F(
			TString::Format("hqspta%d", i),
			TString::Format("Q with single PPAC(%d) track reaction point", i),
			150, -10, 5
		);
	}
	// output tree
	TTree opt("tree", "Q in different condition");
	// output data
	// multiple PPAC tracking target x and y
	double mptx, mpty;
	// single PPAC tracking approximate target x and y
	double sptax[3], sptay[3];
	// flag of single PPAC tracking
	int spt_flag;
	// q value from multiple PPAC tracking
	double q_ppac;
	// q value from single PPAC tracking
	double q_spta[3];
	// setup output branches
	opt.Branch("mptx", &mptx, "mptx/D");
	opt.Branch("mpty", &mpty, "mpty/D");
	opt.Branch("spt_flag", &spt_flag, "spt_flag/I");
	opt.Branch("sptax", sptax, "sptax[3]/D");
	opt.Branch("sptay", sptay, "sptay[3]/D");
	opt.Branch("q_mpt", &q_ppac, "qmpt/D");
	opt.Branch("q_spta", q_spta, "qspta[3]/D");

	double c14_kinetic, d_kinetic, c15_kinetic;

	// total number of entries
	long long entries = ipt->GetEntries();
	for (long long entry = 0; entry < entries; ++entry) {
		// get data
		ipt->GetEntry(entry);

		if (event.ppac_xflag <= 2 || event.ppac_xflag == 4) continue;
		if (event.ppac_yflag <= 2 || event.ppac_yflag == 4) continue;

		// calculate reaction point from PPAC data
		// initialize correct parameters
		double ppac_correct[2][3] = {
			{0.0, 0.0, 0.0},
			{0.0, 0.0, 0.0}
		};
		// calculate correct offsets
		PpacOffsetX(ppac_correct[0][0], ppac_correct[0][1], ppac_correct[0][2]);
		PpacOffsetY(ppac_correct[1][0], ppac_correct[1][1], ppac_correct[1][2]);
		// calculate reaction point
		double ppac_cx[3], ppac_cy[3];
		// correct PPAC
		for (int i = 0; i < 3; ++i) {
			ppac_cx[i] = event.ppac_x[i] - ppac_correct[0][i];
			ppac_cy[i] = event.ppac_y[i] - ppac_correct[1][i];
		}
		// slope and intercept
		double xk, yk;
		// fit and get reaction point
		TrackMultiplePpac(event.ppac_xflag, ppac_xz, ppac_cx, xk, mptx);
		TrackMultiplePpac(event.ppac_yflag, ppac_yz, ppac_cy, yk, mpty);
		q_ppac = PdProcess(
			event, mptx, mpty,
			c14_kinetic, d_kinetic, c15_kinetic
		);
		hist_q_mpt.Fill(q_ppac);

		// single PPAC tracking
		spt_flag = 0;
		for (int i = 0; i < 3; ++i) {
			unsigned short flag = 1 << i;
			if ((event.ppac_xflag & event.ppac_yflag & flag) != flag) {
				continue;
			}
			spt_flag |= flag;
			// apprximate method
			sptax[i] = ApproximateTrack(
				c14_kinetic, d_kinetic, ppac_xz[i],
				event.c_x[0], event.d_x, event.d_y, ppac_cx[i]
			);
			sptay[i] = ApproximateTrack(
				c14_kinetic, d_kinetic, ppac_yz[i],
				event.c_y[0], event.d_y, event.d_x, ppac_cy[i]
			);
			q_spta[i] = PdProcess(
				event, sptax[i], sptay[i],
				c14_kinetic, d_kinetic, c15_kinetic
			);
			hist_q_spta[i]->Fill(q_spta[i]);

			hist_dx[i]->Fill(sptax[i] - mptx);
			hist_dy[i]->Fill(sptay[i] - mpty);
		}

		// fill
		opt.Fill();
	}

	// save histograms
	hist_q_mpt.Write();
	for (int i = 0; i < 3; ++i) hist_q_spta[i]->Write();
	for (int i = 0; i < 3; ++i) hist_dx[i]->Write();
	for (int i = 0; i < 3; ++i) hist_dy[i]->Write();
	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}