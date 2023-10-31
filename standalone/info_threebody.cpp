#include <iostream>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/channel_event.h"
#include "include/event/dssd_event.h"
#include "include/event/particle_event.h"
#include "include/event/particle_type_event.h"
#include "include/event/ta_event.h"
#include "include/event/t0_event.h"

using namespace ribll;

int main() {
	// output file name
	TString output_file_name = TString::Format(
		"%sinfo/threebody.root",
		kGenerateDataPath
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "information of three body");
	// output data
	// indexes and layers
	int csi_index;
	int layer[2];
	// energy
	double tafd_energy;
	// double csi_energy;
	double t0_energy[2];
	// double d1_energy[2];
	// double d2_energy[2];
	// double d3_energy[2];
	// double ssd_energy[3];
	double csi_channel;
	double d1_channel[2];
	double d2_channel[2];
	double d3_channel[2];
	double ssd_channel[3];
	// double d1_x_channel[2];
	// double d1_y_channel[2];
	// double d2_x_channel[2];
	// double d2_y_channel[2];
	// double d3_x_channel[2];
	// double d3_y_channel[2];
	// time
	// double d1_x_time[2];
	// double d1_y_time[2];
	// double d2_x_time[2];
	// double d2_y_time[2];
	// double d3_x_time[2];
	// double d3_y_time[2];
	// double r_time;
	// position
	double d1_x[2];
	double d1_y[2];
	double d2_x[2];
	double d2_y[2];
	double d3_x[2];
	double d3_y[2];
	// recoil position
	double rx;
	double ry;
	// target position
	double tx;
	double ty;
	// strips
	// unsigned int d1_x_strip[2];
	// unsigned int d1_y_strip[2];
	// unsigned int d2_x_strip[2];
	// unsigned int d2_y_strip[2];
	// unsigned int d3_x_strip[2];
	// unsigned int d3_y_strip[2];
	// unsigned int r_ring_strip[2];
	// unsigned int r_pi_strip[2];
	// state
	// int state;
	double out_q;

	// setup output branches
	// indexes and layers
	opt.Branch("csi_index", &csi_index, "ci/I");
	opt.Branch("layer", layer, "layer[2]/I");
	// energy
	opt.Branch("tafd_energy", &tafd_energy, "tafde/D");
	// opt.Branch("csi_energy", &csi_energy, "csie/D");
	opt.Branch("t0_energy", t0_energy, "t0e[2]/D");
	// opt.Branch("d1_energy", d1_energy, "d1e[2]/D");
	// opt.Branch("d2_energy", d2_energy, "d2e[2]/D");
	// opt.Branch("d3_energy", d3_energy, "d3e[2]/D");
	// opt.Branch("ssd_energy", ssd_energy, "ssde[3]/D");
	opt.Branch("csi_channel", &csi_channel, "csic/D");
	opt.Branch("d1_channel", d1_channel, "d1c[2]/D");
	opt.Branch("d2_channel", d2_channel, "d2c[2]/D");
	opt.Branch("d3_channel", d3_channel, "d3c[2]/D");
	opt.Branch("ssd_channel", ssd_channel, "ssdc[3]/D");
	// opt.Branch("d1_x_channel", d1_x_channel, "d1xc[2]/D");
	// opt.Branch("d1_y_channel", d1_y_channel, "d1yc[2]/D");
	// opt.Branch("d2_x_channel", d2_x_channel, "d2xc[2]/D");
	// opt.Branch("d2_y_channel", d2_y_channel, "d2yc[2]/D");
	// opt.Branch("d3_x_channel", d3_x_channel, "d3xc[2]/D");
	// opt.Branch("d3_y_channel", d3_y_channel, "d3yc[2]/D");
	// time
	// opt.Branch("d1_x_time", d1_x_time, "d1xt[2]/D");
	// opt.Branch("d1_y_time", d1_y_time, "d1yt[2]/D");
	// opt.Branch("d2_x_time", d2_x_time, "d2xt[2]/D");
	// opt.Branch("d2_y_time", d2_y_time, "d2yt[2]/D");
	// opt.Branch("d3_x_time", d3_x_time, "d3xt[2]/D");
	// opt.Branch("d3_y_time", d3_y_time, "d3yt[2]/D");
	// opt.Branch("r_time", &r_time, "rt/D");
	// position
	opt.Branch("d1x", d1_x, "d1x[2]/D");
	opt.Branch("d1y", d1_y, "d1y[2]/D");
	opt.Branch("d2x", d2_x, "d2x[2]/D");
	opt.Branch("d2y", d2_y, "d2y[2]/D");
	opt.Branch("d3x", d3_x, "d3x[2]/D");
	opt.Branch("d3y", d3_y, "d3y[2]/D");
	opt.Branch("rx", &rx, "rx/D");
	opt.Branch("ry", &ry, "ry/D");
	opt.Branch("tx", &tx, "tx/D");
	opt.Branch("ty", &ty, "ty/D");
	// opt.Branch("d1xs", d1_x_strip, "d1xs[2]/D");
	// opt.Branch("d1ys", d1_y_strip, "d1ys[2]/D");
	// opt.Branch("d2xs", d2_x_strip, "d2xs[2]/D");
	// opt.Branch("d2ys", d2_y_strip, "d2ys[2]/D");
	// opt.Branch("rrs", &r_ring_strip, "rrs/D");
	// opt.Branch("rps", &r_pi_strip, "rps/D");
	// opt.Branch("state", &state, "state/I");
	opt.Branch("q", &out_q, "q/D");

	for (unsigned int run = 618; run <= 652; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		std::cout << "Processing run " << run << "\n";
		// input file name
		TString channel_file_name = TString::Format(
			"%s%sC14-10Be-4He-2H-%04d.root",
			kGenerateDataPath,
			kChannelDir,
			run
		);
		// input file
		TFile channel_file(channel_file_name, "read");
		// input tree
		TTree *ipt = (TTree*)channel_file.Get("tree");
		if (!ipt) {
			std::cerr << "Error: Get tree from "
				<< channel_file_name << " failed.\n";
			return -1;
		}
		// input channel event
		ChannelEvent channel;
		// t0 index
		int t0_index[8];
		// setup input branches
		channel.SetupInput(ipt);
		ipt->SetBranchAddress("t0_index", t0_index);

		std::vector<long long> valid_entries;
		std::vector<int> taf_indexes;
		std::vector<int> be10_indexes;
		std::vector<int> he4_indexes;
		std::vector<double> t0_energies[2];
		std::vector<double> q_values;

		// loop to record information in channel
		for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
			ipt->GetEntry(entry);
			valid_entries.push_back(channel.entry);
			taf_indexes.push_back(channel.taf_index);
			be10_indexes.push_back(t0_index[0]);
			he4_indexes.push_back(t0_index[1]);
			t0_energies[0].push_back(channel.daughter_energy[0]);
			t0_energies[1].push_back(channel.daughter_energy[1]);
			double q = channel.daughter_energy[0] + channel.daughter_energy[1]
				+ channel.recoil_energy - channel.parent_energy;
			q_values.push_back(q);
			// double q = channel.parent_energy + 1.006
			// 	- channel.daughter_energy[0] - channel.daughter_energy[1];
			// if (q > -3.0 && q < 2.5) states.push_back(0);
			// else if (q > 3.0 && q < 9.0) states.push_back(1);
			// else states.push_back(-1);
		}
		// close file
		channel_file.Close();

		// open telescope file to get more information
		TString t0_tele_file_name = TString::Format(
			"%s%st0-telescope-ta-%04u.root",
			kGenerateDataPath,
			kTelescopeDir,
			run
		);
		// t0 file
		TFile t0_tele_file(t0_tele_file_name, "read");
		// t0 tree
		TTree *t0_tree = (TTree*)t0_tele_file.Get("tree");
		if (!t0_tree) {
			std::cerr << "Error: Get tree from "
				<< t0_tele_file_name << " failed.\n";
			return -1;
		}
		// add T0 particle identify file
		t0_tree->AddFriend("pid=tree", TString::Format(
			"%s%st0-particle-type-ta-%04u.root",
			kGenerateDataPath,
			kParticleIdentifyDir,
			run
		));
		// // add T0D1 normalize result friend
		// t0_tree->AddFriend("t0d1=tree", TString::Format(
		// 	"%s%st0d1-result-ta-%04u.root",
		// 	kGenerateDataPath,
		// 	kNormalizeDir,
		// 	run
		// ));
		// // add T0D2 normalize result friend
		// t0_tree->AddFriend("t0d2=tree", TString::Format(
		// 	"%s%st0d2-result-ta-%04u.root",
		// 	kGenerateDataPath,
		// 	kNormalizeDir,
		// 	run
		// ));
		for (int i = 0; i < 6; ++i) {
			// add TAF telescope tree as friend
			t0_tree->AddFriend(
				TString::Format("taf%d=tree", i),
				TString::Format(
					"%s%staf%d-telescope-ta-%04u.root",
					kGenerateDataPath,
					kTelescopeDir,
					i,
					run
				)
			);
			// add TAF ADSSD fundamental tree as friend
			// t0_tree->AddFriend(
			// 	TString::Format("tafd%d=tree", i),
			// 	TString::Format(
			// 		"%s%stafd%d-fundamental-ta-%04u.root",
			// 		kGenerateDataPath,
			// 		kFundamentalDir,
			// 		i,
			// 		run
			// 	)
			// );
		}
		// add xppac as friend
		t0_tree->AddFriend(
			"xppac=tree",
			TString::Format(
				"%s%sxppac-particle-ta-%04u.root",
				kGenerateDataPath,
				kParticleDir,
				run
			)
		);
		// input events
		// T0 telescope event
		T0Event t0;
		ParticleTypeEvent t0_type;
		// T0D1 normalize result event
		// DssdFundamentalEvent t0d1;
		// T0D2 normalize result event
		// DssdFundamentalEvent t0d2;
		// TAF telescope events
		TaEvent taf[6];
		// TAF ADSSD fundamental events
		// DssdFundamentalEvent tafd[6];
		// XPPAC events
		ParticleEvent xppac;

		// setup input branches
		t0.SetupInput(t0_tree);
		t0_type.SetupInput(t0_tree, "pid.");
		// t0d1.SetupInput(t0_tree, "t0d1.");
		// t0d2.SetupInput(t0_tree, "t0d2.");
		for (int i = 0; i < 6; ++i) {
			taf[i].SetupInput(t0_tree, TString::Format("taf%d.", i).Data());
			// tafd[i].SetupInput(t0_tree, TString::Format("tafd%d.", i).Data());
		}
		xppac.SetupInput(t0_tree, "xppac.");

		for (size_t i = 0; i < valid_entries.size(); ++i) {
			t0_tree->GetEntry(valid_entries[i]);
			// energy
			tafd_energy = taf[taf_indexes[i]].energy[0][0];
			t0_energy[0] = t0_energies[0][i];
			t0_energy[1] = t0_energies[1][i];
			if (taf[taf_indexes[i]].flag[0] == 0x3) {
				csi_index = taf_indexes[i]*2;
			} else if (taf[taf_indexes[i]].flag[0] == 0x5) {
				csi_index = taf_indexes[i]*2 + 1;
			} else {
				csi_index = -1;
			}
			csi_channel = taf[taf_indexes[i]].energy[0][1];
			layer[0] = t0_type.layer[be10_indexes[i]];
			layer[1] = t0_type.layer[he4_indexes[i]];
			d1_channel[0] = t0.energy[be10_indexes[i]][0];
			d1_channel[1] = t0.energy[he4_indexes[i]][0];
			d2_channel[0] = t0.energy[be10_indexes[i]][1];
			d2_channel[1] = t0.energy[he4_indexes[i]][1];
			d3_channel[0] = t0.energy[be10_indexes[i]][2];
			d3_channel[1] = t0.energy[he4_indexes[i]][2];
			ssd_channel[0] = t0.ssd_energy[0];
			ssd_channel[1] = t0.ssd_energy[1];
			ssd_channel[2] = t0.ssd_energy[2];
			// position
			d1_x[0] = t0.x[be10_indexes[i]][0];
			d1_y[0] = t0.y[be10_indexes[i]][0];
			d2_x[0] = t0.x[be10_indexes[i]][1];
			d2_y[0] = t0.y[be10_indexes[i]][1];
			d3_x[0] = t0.x[be10_indexes[i]][2];
			d3_y[0] = t0.y[be10_indexes[i]][2];
			d1_x[1] = t0.x[he4_indexes[i]][0];
			d1_y[1] = t0.y[he4_indexes[i]][0];
			d2_x[1] = t0.x[he4_indexes[i]][1];
			d2_y[1] = t0.y[he4_indexes[i]][1];
			d3_x[1] = t0.x[he4_indexes[i]][2];
			d3_y[1] = t0.y[he4_indexes[i]][2];
			ROOT::Math::Polar3DVector recoil_position(
				taf[taf_indexes[i]].radius[0],
				taf[taf_indexes[i]].theta[0],
				taf[taf_indexes[i]].phi[0]
			);
			rx = recoil_position.X();
			ry = recoil_position.Y();
			tx = xppac.x[3];
			ty = xppac.y[3];
			out_q = q_values[i];

			opt.Fill();
		}
		t0_tele_file.Close();
	}
	opf.cd();
	// save trees
	opt.Write();
	// close files
	opf.Close();
	return 0;
}