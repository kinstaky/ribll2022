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
#include "include/event/threebody_info_event.h"

using namespace ribll;

int FillResult(
	unsigned short flag,
	int result_hit,
	double *result_channel,
	double *result_time,
	unsigned short *result_strip,
	int &hit,
	double *channel,
	double *time,
	unsigned int *strip
) {
	hit = 0;
	for (int i = 0; i < 8; ++i) {
		if ((flag & (1 << i)) == 0) continue;
		if (i >= result_hit) {
			std::cerr << "Error: Found flag bit over result hit.\n";
			return -1;
		}
		channel[hit] = result_channel[i];
		time[hit] = result_time[i];
		strip[hit] = result_strip[i];
		++hit;
	}
	return 0;
}

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
	ThreeBodyInfoEvent event;
	// setup output branches
	event.SetupOutput(&opt);

	for (unsigned int run = 618; run <= 716; ++run) {
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
		// identify flag
		std::vector<bool> *identify =
			(std::vector<bool>*)t0_tele_file.Get("identify");
		if (!identify || !(identify->at(0))) {
			// add T0 particle identify file
			t0_tree->AddFriend("pid=tree", TString::Format(
				"%s%st0-particle-type-ta-%04u.root",
				kGenerateDataPath,
				kParticleIdentifyDir,
				run
			));
		}
		// add T0D1 normalize result friend
		t0_tree->AddFriend("t0d1=tree", TString::Format(
			"%s%st0d1-result-ta-%04u.root",
			kGenerateDataPath,
			kNormalizeDir,
			run
		));
		// add T0D2 normalize result friend
		t0_tree->AddFriend("t0d2=tree", TString::Format(
			"%s%st0d2-result-ta-%04u.root",
			kGenerateDataPath,
			kNormalizeDir,
			run
		));
		// add T0D3 normalize result friend
		t0_tree->AddFriend("t0d3=tree", TString::Format(
			"%s%st0d3-result-ta-%04u.root",
			kGenerateDataPath,
			kNormalizeDir,
			run
		));
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
			t0_tree->AddFriend(
				TString::Format("tafd%d=tree", i),
				TString::Format(
					"%s%stafd%d-fundamental-ta-%04u.root",
					kGenerateDataPath,
					kFundamentalDir,
					i,
					run
				)
			);
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
		// T0Dx normalize result event
		DssdFundamentalEvent dssd_result[3];
		// TAF telescope events
		TaEvent taf[6];
		// TAF ADSSD fundamental events
		DssdFundamentalEvent tafd[6];
		// XPPAC events
		ParticleEvent xppac;

		// setup input branches
		t0.SetupInput(t0_tree);
		if (!identify || !(identify->at(0))) {
			t0_type.SetupInput(t0_tree, "pid.");
		} else {
			t0_tree->SetBranchAddress("mass", t0_type.mass);
			t0_tree->SetBranchAddress("charge", t0_type.charge);
			t0_tree->SetBranchAddress("layer", t0_type.layer);
		}
		dssd_result[0].SetupInput(t0_tree, "t0d1.");
		dssd_result[1].SetupInput(t0_tree, "t0d2.");
		dssd_result[2].SetupInput(t0_tree, "t0d3.");
		for (int i = 0; i < 6; ++i) {
			taf[i].SetupInput(t0_tree, TString::Format("taf%d.", i).Data());
			tafd[i].SetupInput(t0_tree, TString::Format("tafd%d.", i).Data());
		}
		xppac.SetupInput(t0_tree, "xppac.");
		t0_tree->SetBranchAddress("xppac.xflag", &event.ppac_xflag);
		t0_tree->SetBranchAddress("xppac.yflag", &event.ppac_yflag);

		for (size_t i = 0; i < valid_entries.size(); ++i) {
			t0_tree->GetEntry(valid_entries[i]);
			// energy
			event.tafd_energy = taf[taf_indexes[i]].energy[0][0];
			event.t0_energy[0] = t0_energies[0][i];
			event.t0_energy[1] = t0_energies[1][i];
			if (taf[taf_indexes[i]].flag[0] == 0x3) {
				event.csi_index = taf_indexes[i]*2;
			} else if (taf[taf_indexes[i]].flag[0] == 0x5) {
				event.csi_index = taf_indexes[i]*2 + 1;
			} else {
				event.csi_index = -1;
			}
			event.csi_channel = taf[taf_indexes[i]].energy[0][1];
			event.layer[0] = t0_type.layer[be10_indexes[i]];
			event.layer[1] = t0_type.layer[he4_indexes[i]];

			for (int j = 0; j < 3; ++j) {
				event.be_channel[j] = t0.energy[be10_indexes[i]][j];
				event.he_channel[j] = t0.energy[he4_indexes[i]][j];
			}
			for (int j = 0; j < 3; ++j) {
				event.ssd_channel[j] = t0.ssd_energy[j];
			}
			// position
			for (int j = 0; j < 3; ++j) {
				event.be_x[j] = t0.x[be10_indexes[i]][j];
				event.be_y[j] = t0.y[be10_indexes[i]][j];
				event.he_x[j] = t0.x[he4_indexes[i]][j];
				event.he_y[j] = t0.y[he4_indexes[i]][j];
			}
			ROOT::Math::Polar3DVector recoil_position(
				taf[taf_indexes[i]].radius[0],
				taf[taf_indexes[i]].theta[0],
				taf[taf_indexes[i]].phi[0]
			);
			event.d_x = recoil_position.X();
			event.d_y = recoil_position.Y();
			event.tx = xppac.x[3];
			event.ty = xppac.y[3];
			for (int j = 0; j < 3; ++j) {
				event.ppac_x[j] = xppac.x[j];
				event.ppac_y[j] = xppac.y[j];
			}


			// fill t0d1 10Be result
			if (FillResult(
				t0.dssd_flag[be10_indexes[i]][0] >> 8,
				dssd_result[0].back_hit,
				dssd_result[0].back_energy,
				dssd_result[0].back_time,
				dssd_result[0].back_strip,
				event.be_x_hit[0],
				event.be_x_channel[0],
				event.be_x_time[0],
				event.be_x_strip[0]
			)) {
				std::cerr << "Error: Run " << run
					<< ", entry " << valid_entries[i]
					<< " 10Be T0D1 back. \n";
				return -1;
			}
			if (FillResult(
				t0.dssd_flag[be10_indexes[i]][0],
				dssd_result[0].front_hit,
				dssd_result[0].front_energy,
				dssd_result[0].front_time,
				dssd_result[0].front_strip,
				event.be_y_hit[0],
				event.be_y_channel[0],
				event.be_y_time[0],
				event.be_y_strip[0]
			)) {
				std::cerr << "Error: Run " << run
					<< ", entry " << valid_entries[i]
					<< " 10Be T0D1 front. \n";
				return -1;
			}

			// fill t0d1 4He result
			if (FillResult(
				t0.dssd_flag[he4_indexes[i]][0] >> 8,
				dssd_result[0].back_hit,
				dssd_result[0].back_energy,
				dssd_result[0].back_time,
				dssd_result[0].back_strip,
				event.he_x_hit[0],
				event.he_x_channel[0],
				event.he_x_time[0],
				event.he_x_strip[0]
			)) {
				std::cerr << "Error: Run " << run
					<< ", entry " << valid_entries[i]
					<< " 4He T0D1 back. \n";
				return -1;
			}
			if(FillResult(
				t0.dssd_flag[he4_indexes[i]][0],
				dssd_result[0].front_hit,
				dssd_result[0].front_energy,
				dssd_result[0].front_time,
				dssd_result[0].front_strip,
				event.he_y_hit[0],
				event.he_y_channel[0],
				event.he_y_time[0],
				event.he_y_strip[0]
			)) {
				std::cerr << "Error: Run " << run
					<< ", entry " << valid_entries[i]
					<< " 4He T0D1 front. \n";
				return -1;
			}

			// fill T0D2 and T0D3 10Be result
			for (int j = 1; j < event.layer[0]+1; ++j) {
				if (j >= 3) break;
				if (FillResult(
					t0.dssd_flag[be10_indexes[i]][j],
					dssd_result[j].front_hit,
					dssd_result[j].front_energy,
					dssd_result[j].front_time,
					dssd_result[j].front_strip,
					event.be_x_hit[j],
					event.be_x_channel[j],
					event.be_x_time[j],
					event.be_x_strip[j]
				)) {
					std::cerr << "Error: Run " << run
						<< ", entry " << valid_entries[i]
						<< " 10Be T0D" << j+1 << " front. \n";
					return -1;
				}
				if (FillResult(
					t0.dssd_flag[be10_indexes[i]][j] >> 8,
					dssd_result[j].back_hit,
					dssd_result[j].back_energy,
					dssd_result[j].back_time,
					dssd_result[j].back_strip,
					event.be_y_hit[j],
					event.be_y_channel[j],
					event.be_y_time[j],
					event.be_y_strip[j]
				)) {
					std::cerr << "Error: Run " << run
						<< ", entry " << valid_entries[i]
						<< " 10Be T0D" << j+1 << " back. \n";
					return -1;
				}
			}

			// fill T0D2 and T0D3 4He result
			for (int j = 1; j < event.layer[1]+1; ++j) {
				if (j >= 3) break;
				// 4He result
				if (FillResult(
					t0.dssd_flag[he4_indexes[i]][j],
					dssd_result[j].front_hit,
					dssd_result[j].front_energy,
					dssd_result[j].front_time,
					dssd_result[j].front_strip,
					event.he_x_hit[j],
					event.he_x_channel[j],
					event.he_x_time[j],
					event.he_x_strip[j]
				)) {
					std::cerr << "Error: Run " << run
						<< ", entry " << valid_entries[i]
						<< " 4He T0D" << j+1 << " front. \n";
					return -1;
				}
				if (FillResult(
					t0.dssd_flag[he4_indexes[i]][j] >> 8,
					dssd_result[j].back_hit,
					dssd_result[j].back_energy,
					dssd_result[j].back_time,
					dssd_result[j].back_strip,
					event.he_y_hit[j],
					event.he_y_channel[j],
					event.he_y_time[j],
					event.he_y_strip[j]
				)) {
					std::cerr << "Error: Run " << run
						<< ", entry " << valid_entries[i]
						<< " 4He T0D" << j+1 << " back. \n";
					return -1;
				}
			}

			event.d_x_channel = tafd[taf_indexes[i]].front_energy[0];
			event.d_y_channel = tafd[taf_indexes[i]].back_energy[0];
			event.d_x_time = tafd[taf_indexes[i]].front_time[0];
			event.d_y_time = tafd[taf_indexes[i]].back_time[0];
			event.d_x_strip = tafd[taf_indexes[i]].front_strip[0];
			event.d_y_strip = tafd[taf_indexes[i]].back_strip[0];

			event.q = q_values[i];
			event.hole[0] = t0.hole[be10_indexes[i]];
			event.hole[1] = t0.hole[he4_indexes[i]];
			event.run = run;

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