#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TGraph.h>
#include <Math/Vector3D.h>

#include "include/event/threebody_info_event.h"
#include "include/event/channel_v2_event.h"

using namespace ribll;

// Q value correct
constexpr double q_tafd_correct[6] = {
	0.3, -0.3, -0.4, 0.0, 0.4, 0.4
};
constexpr double q_csi_correct[12] = {
	-0.46, -0.25, 0.18, 0.49,
	0.49, 0.55, 0.28, 0.30,
	-0.04, -0.10, -0.12, -0.64
};
const double energy_base[3] = {12.02, 15.39, 18.20};


int main() {
	// spectrum V3 file
	TString input_file_name = TString::Format(
		"%s%sC14-10Be-4He-2H-v2-sim.root", kGenerateDataPath, kChannelDir
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
		"%s%sefficiency-v2.root",
		kGenerateDataPath,
		kSimulateDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// total efficiency
	TGraph graph_efficiency[3];
	TTree opt("tree", "efficiency");
	int valid;
	opt.Branch("valid", &valid, "v/I");

	int valid_count[3][100];
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 100; ++j) {
			valid_count[i][j] = 0;
		}
	}

	// loop to process
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		// get data
		ipt->GetEntry(entry);

		if (channel.hole != 0) continue;
		// get flags
		valid = 0;
		if (!channel.t0_valid) valid |= 1;
		if (!channel.tafcsi_valid && channel.tafd_front_strip >= 13) {
		} else if (channel.tafcsi_valid && channel.tafd_front_strip <= 13) {
		} else {
			valid |= 2;
		}

		if (!channel.ppac_valid) valid |= 4;
		if (pow(channel.tx, 2.0) + pow(channel.ty, 2.0) >= 196.0) {
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
		double be_kinetic = channel.fragment_kinetic[0];
		// 4He kinetic energy
		double he_kinetic = channel.fragment_kinetic[1];
		// deutron kinetic energy
		double d_kinetic = channel.recoil_kinetic;

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
		double c_momentum = p_c.R();
		// 14C kinematic energy
		double c_kinetic =
			sqrt(pow(c_momentum, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

		// three body Q value
		double q = be_kinetic + he_kinetic + d_kinetic - c_kinetic;

		int be_state = -1;
		// get 10Be state from Q value
		if (q < -11.0 && q > -13.0) be_state = 0;
		else if (q < -14.5 && q > -16.3) be_state = 1;
		else if (q < -17.0 && q > -20.0) be_state = 2;
		// else if (q < -19 && q > -20.5) be_state = 3;
		else be_state = -1;


		if (c_kinetic < 360.0) valid |= 16;

		// index
		int i1 = channel.entry / 10'000'000;
		int i2 = (channel.entry % 10'000'000) / 100000;
		if (i1 == be_state && valid == 0) {
			++valid_count[i1][i2];
		}
		opt.Fill();
	}

	// fill to graph
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 100; ++j) {
			graph_efficiency[i].AddPoint(
				energy_base[i] + 0.2 * j,
				valid_count[i][j] / 100000.0
			);
		}
	}

	for (int i = 0; i < 3; ++i) {
		graph_efficiency[i].Write(TString::Format("g%d", i));
	}
	opt.Write();
	opf.Close();
	ipf.Close();
	return 0;
}