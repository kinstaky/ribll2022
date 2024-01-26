#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/threebody_info_event.h"
#include "include/ppac_track.h"

using namespace ribll;


double ThreeBodyProcess(
	const ThreeBodyInfoEvent &event,
	const double ppac_parameter_x,
	const double ppac_parameter_y,
	double &be_kinetic,
	double &he_kinetic,
	double &d_kinetic,
	double &c_kinetic
) {
	// 10Be kinetic energy
	// T0D1 energy
	be_kinetic = t0_param[0][0] + t0_param[0][1] * event.be_channel[0];
	// T0D2 energy
	be_kinetic += t0_param[1][0] + t0_param[1][1] * event.be_channel[1];
	// T0D3 energy
	if (event.layer[0] > 1) {
		be_kinetic += t0_param[2][0] + t0_param[2][1] * event.be_channel[2];
	}
	// T0S1 energy
	if (event.layer[0] > 2) {
		be_kinetic += t0_param[3][0] + t0_param[3][1] * event.ssd_channel[0];
	}
	// T0S2 energy
	if (event.layer[0] > 3) {
		be_kinetic += t0_param[4][0] + t0_param[4][1] * event.ssd_channel[1];
	}
	// T0S3 energy
	if (event.layer[0] > 4) {
		be_kinetic += t0_param[5][0] + t0_param[5][1] * event.ssd_channel[2];
	}

	// 4He kinetic energy
	// T0D1 energy
	he_kinetic = t0_param[0][0] + t0_param[0][1] * event.he_channel[0];
	// T0D2 energy
	he_kinetic += t0_param[1][0] + t0_param[1][1] * event.he_channel[1];
	// T0D3 energy
	if (event.layer[1] > 1) {
		he_kinetic += t0_param[2][0] + t0_param[2][1] * event.he_channel[2];
	}
	// T0S1 energy
	if (event.layer[1] > 2) {
		he_kinetic += t0_param[3][0] + t0_param[3][1] * event.ssd_channel[0];
	}
	// T0S2 energy
	if (event.layer[1] > 3) {
		he_kinetic += t0_param[4][0] + t0_param[4][1] * event.ssd_channel[1];
	}
	// T0S3 energy
	if (event.layer[1] > 4) {
		he_kinetic += t0_param[5][0] + t0_param[5][1] * event.ssd_channel[2];
	}

	double a0 = csi_param[event.csi_index][0];
	double a1 = csi_param[event.csi_index][1];
	double a2 = csi_param[event.csi_index][2];
	// 2H kinetic energy
	d_kinetic = event.tafd_energy + pow(
		(event.csi_channel - a2) / a0,
		1.0 / a1
	);

	double ppac_correct[2][3] = {
		{ppac_parameter_x, 0.0, 0.0},
		{ppac_parameter_y, 0.0, 0.0}
	};
	PpacOffsetX(ppac_correct[0][0], ppac_correct[0][1], ppac_correct[0][2]);
	PpacOffsetY(ppac_correct[1][0], ppac_correct[1][1], ppac_correct[1][2]);
	// calculate reaction point
	double ppac_cx[3], ppac_cy[3];
	for (int i = 0; i < 3; ++i) {
		ppac_cx[i] = event.ppac_x[i] - ppac_correct[0][i];
		ppac_cy[i] = event.ppac_y[i] - ppac_correct[1][i];
	}
	// slope and intercept
	double xk, yk, xb, yb;
	TrackMultiplePpac(event.ppac_xflag, ppac_xz, ppac_cx, xk, xb);
	TrackMultiplePpac(event.ppac_yflag, ppac_yz, ppac_cy, yk, yb);

	// 10Be momentum
	double be_momentum = MomentumFromKinetic(mass_10be, be_kinetic);
	// 10Be momentum vector
	ROOT::Math::XYZVector p_be(
		event.be_x[0] - xb,
		event.be_y[0] - yb,
		100.0
	);
	p_be = p_be.Unit() * be_momentum;

	// 4He momentum
	double he_momentum = MomentumFromKinetic(mass_4he, he_kinetic);
	// 4He momentum vector
	ROOT::Math::XYZVector p_he(
		event.he_x[0] - xb,
		event.he_y[0] - yb,
		100.0
	);
	p_he = p_he.Unit() * he_momentum;

	// 2H momentum
	double d_momentum = MomentumFromKinetic(mass_2h, d_kinetic);
	// 2H momentum vector
	ROOT::Math::XYZVector p_d(
		event.d_x - xb,
		event.d_y - yb,
		135.0
	);
	p_d = p_d.Unit() * d_momentum;

	// beam 14C momentum vector
	ROOT::Math::XYZVector p_c = p_be + p_he + p_d;
// std::cout << p_be.X() << ", " << p_be.Y() << ", " << p_be.Z() << "\n"
// 	<< p_he.X() << ", " << p_he.Y() << ", " << p_he.Z() << "\n"
// 	<< p_d.X() << ", " << p_d.Y() << ", " << p_d.Z() << "\n"
// 	<< p_c.X() << ", " << p_c.Y() << ", " << p_c.Z() << "\n";

	// 14C momentum
	double c_momentum = p_c.R();
	// 14C kinetic energy
	c_kinetic =
		sqrt(pow(c_momentum, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

	// three-fold Q value
	double q = be_kinetic + he_kinetic + d_kinetic - c_kinetic;

	return q;
}


int main() {
	// input file name
	TString input_file_name = TString::Format(
		"%s%sthreebody.root", kGenerateDataPath, kInformationDir
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
	// add CsI group result
	ipt->AddFriend("group=tree", TString::Format(
		"%s%sthreebody-group.root", kGenerateDataPath, kOptimizeDir
	));
	// input data
	ThreeBodyInfoEvent event;
	double csi_mean;
	double csi_sigma;
	// setup input branches
	event.SetupInput(ipt);
	ipt->SetBranchAddress("group.q_mean", &csi_mean);
	ipt->SetBranchAddress("group.q_sigma", &csi_sigma);

	constexpr int groups = 64;

	// output file name
	TString output_file_name = TString::Format(
		"%s%sppac-group.root", kGenerateDataPath, kOptimizeDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// group parameters Q value
	std::vector<TH1F> hist_group_q;
	// Q value spectrum for each group
	std::vector<TH1F> hist_q_spectrum;
	for (int i = 0; i < groups; ++i) {
		hist_q_spectrum.emplace_back(
			TString::Format("hq%d", i), "Q",
			100, -23, -8
		);
	}
	// output tree
	TTree opt("tree", "ppac group");
	// output data
	int num;
	double ppac_parameter_x[groups];
	double ppac_parameter_y[groups];
	double group_q[groups];
	double q_mean, q_sigma;
	// setup output branches
	opt.Branch("num", &num, "num/I");
	opt.Branch("lx", ppac_parameter_x, "lx[num]/D");
	opt.Branch("ly", ppac_parameter_y, "ly[num]/D");
	opt.Branch("q", group_q, "q[num]/D");
	opt.Branch("q_mean", &q_mean, "q_mean/D");
	opt.Branch("q_sigma", &q_sigma, "q_sigma/D");

	double q, be_kinetic, he_kinetic, d_kinetic, c_kinetic;

	// total number of entries
	long long entries = ipt->GetEntries();
	for (long long entry = 0; entry < entries; ++entry) {
		// get data
		ipt->GetEntry(entry);

		hist_group_q.emplace_back(
			TString::Format("hgq%lld", entry), "group Q",
			1000, -23, -8
		);

		num = 0;
		if (csi_sigma > 0.2) {
			opt.Fill();
			continue;
		}

		for (double lx = -4.0; lx < 4.0; lx += 1.0) {
			for (double ly = -4.0; ly < 4.0; ly += 1.0) {
				q = ThreeBodyProcess(
					event,
					lx, ly,
					be_kinetic, he_kinetic, d_kinetic, c_kinetic
				);
				// fill histogram
				hist_group_q[entry].Fill(q);
				hist_q_spectrum[num].Fill(q);
				// fill tree
				ppac_parameter_x[num] = lx;
				ppac_parameter_y[num] = ly;
				group_q[num] = q;
				++num;
			}
		}
		q_mean = hist_group_q[entry].GetMean();
		q_sigma = hist_group_q[entry].GetStdDev();
		opt.Fill();
	}

	// save histograms
	for (auto &hist : hist_group_q) hist.Write();
	for (auto &hist : hist_q_spectrum) hist.Write();
	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}