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
#include "include/ppac_track.h"
#include "include/calculator/csi_energy_calculator.h"

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


/// @brief rebuild threebody reaction process
/// @param[inout] event input event
/// @returns Q value
///
double ThreeBodyProcess(ThreeBodyInfoEvent &event) {
	// 10Be momentum
	double be_momentum = MomentumFromKinetic(mass_10be, event.t0_energy[0]);
	// 10Be momentum vector
	ROOT::Math::XYZVector p_be(
		event.be_x[0] - event.tx,
		event.be_y[0] - event.ty,
		100.0
	);
	p_be = p_be.Unit() * be_momentum;

	// 4He momentum
	double he_momentum = MomentumFromKinetic(mass_4he, event.t0_energy[1]);
	// 4He momentum vector
	ROOT::Math::XYZVector p_he(
		event.he_x[0] - event.tx,
		event.he_y[0] - event.ty,
		100.0
	);
	p_he = p_he.Unit() * he_momentum;

	// 2H momentum
	double d_momentum = MomentumFromKinetic(mass_2h, event.taf_energy);
	// 2H momentum vector
	ROOT::Math::XYZVector p_d(
		event.d_x - event.tx,
		event.d_y - event.ty,
		135.0
	);
	p_d = p_d.Unit() * d_momentum;

	// beam 14C momentum vector
	ROOT::Math::XYZVector p_c = p_be + p_he + p_d;

	// 14C momentum
	double c_momentum = p_c.R();
	// 14C kinematic energy
	event.c14_kinetic =
		sqrt(pow(c_momentum, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

	// three-fold Q value
	double q = event.t0_energy[0] + event.t0_energy[1]
		+ event.taf_energy - event.c14_kinetic;

	return q;
}


int main(int argc, char **argv) {
	// recoil particle mass number
	int recoil_mass = 2;
	if (argc > 1) {
		recoil_mass = atoi(argv[1]);
	}
	if (recoil_mass <= 0 || recoil_mass > 3) {
		std::cerr << "Error: Invalid recoil particle mass " << recoil_mass << "\n";
		return -1;
	}

	// output file name
	TString output_file_name = TString::Format(
		"%sinfo/threebody%s.root",
		kGenerateDataPath,
		argc == 1 ? "" : TString::Format("-%dH", recoil_mass).Data()
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "information of three body");
	// output data
	ThreeBodyInfoEvent event;
	// int vppac_num, xppac_num;
	int xppac_num;
	// setup output branches
	event.SetupOutput(&opt);
	// opt.Branch("vnum", &vppac_num, "vnum/I");
	opt.Branch("xnum", &xppac_num, "xnum/I");

	// possible 2H stopped in ADSSD
	int possible_2H_num[4];
	for (int i = 0; i < 4; ++i) {
		possible_2H_num[i] = 0;
	}
	// total number of possible 2H events
	int total_possible_2H = 0;

	elc::CsiEnergyCalculator calculator("2H");

	for (unsigned int run = 618; run <= 716; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		std::cout << "Processing run " << run << "\n";
		// input file name
		TString channel_file_name = TString::Format(
			"%s%sC14-10Be-4He-%dH-%04d.root",
			kGenerateDataPath,
			kChannelDir,
			recoil_mass,
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
		// PPAC flag
		int ppac_flag;
		// setup input branches
		channel.SetupInput(ipt);
		ipt->SetBranchAddress("t0_index", t0_index);
		ipt->SetBranchAddress("ppac_flag", &ppac_flag);

		std::vector<long long> valid_entries;
		std::vector<int> taf_indexes;
		std::vector<int> be10_indexes;
		std::vector<int> he4_indexes;
		std::vector<int> ppac_flags;

		// loop to record information in channel
		for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
			ipt->GetEntry(entry);
			valid_entries.push_back(channel.entry);
			taf_indexes.push_back(channel.taf_index);
			be10_indexes.push_back(t0_index[0]);
			he4_indexes.push_back(t0_index[1]);
			ppac_flags.push_back(ppac_flag);
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
		// // add vppac as friend
		// t0_tree->AddFriend(
		// 	"vppac=tree",
		// 	TString::Format(
		// 		"%s%svppac-particle-ta-%04u.root",
		// 		kGenerateDataPath,
		// 		kParticleDir,
		// 		run
		// 	)
		// );
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
		// VPPAC events
		// ParticleEvent vppac;
		// unsigned short vppac_xflag, vppac_yflag;

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
		// vppac.SetupInput(t0_tree, "vppac.");
		// t0_tree->SetBranchAddress("vppac.xflag", &vppac_xflag);
		// t0_tree->SetBranchAddress("vppac.yflag", &vppac_yflag);

		for (size_t i = 0; i < valid_entries.size(); ++i) {
			t0_tree->GetEntry(valid_entries[i]);

			// T0
			// T0 layer
			event.layer[0] = t0_type.layer[be10_indexes[i]];
			event.layer[1] = t0_type.layer[he4_indexes[i]];
			// T0 channels
			for (int j = 0; j < 3; ++j) {
				event.be_channel[j] = t0.energy[be10_indexes[i]][j];
				event.he_channel[j] = t0.energy[he4_indexes[i]][j];
			}
			for (int j = 0; j < 3; ++j) {
				event.ssd_channel[j] = t0.ssd_energy[j];
			}

			// calculate energy
			// 10Be kinetic energy
			// T0D1
			event.t0_energy[0] =
				t0_param[0][0] + t0_param[0][1] * event.be_channel[0];
			// T0D2
			event.t0_energy[0] +=
				t0_param[1][0] + t0_param[1][1] * event.be_channel[1];
			// T0D3
			if (event.layer[0] > 1) {
				event.t0_energy[0] +=
					t0_param[2][0] + t0_param[2][1] * event.be_channel[2];
			}
			// T0S1
			if (event.layer[0] > 2) {
				event.t0_energy[0] +=
					t0_param[3][0] + t0_param[3][1] * event.ssd_channel[0];
			}
			if (event.layer[0] > 3) {
				event.t0_energy[0] +=
					t0_param[4][0] + t0_param[4][1] * event.ssd_channel[1];
			}
			if (event.layer[0] > 4) {
				event.t0_energy[0] +=
					t0_param[5][0] + t0_param[5][1] * event.ssd_channel[2];
			}

			// 4He kinetic energy
			// T0D1
			event.t0_energy[1] =
				t0_param[0][0] + t0_param[0][1] * event.he_channel[0];
			// T0D2
			event.t0_energy[1] +=
				t0_param[1][0] + t0_param[1][1] * event.he_channel[1];
			// T0D3
			if (event.layer[1] > 1) {
				event.t0_energy[1] +=
					t0_param[2][0] + t0_param[2][1] * event.he_channel[2];
			}
			if (event.layer[1] > 2) {
				event.t0_energy[1] +=
					t0_param[3][0] + t0_param[3][1] * event.ssd_channel[0];
			}
			if (event.layer[1] > 3) {
				event.t0_energy[1] +=
					t0_param[4][0] + t0_param[4][1] * event.ssd_channel[1];
			}
			if (event.layer[1] > 4) {
				event.t0_energy[1] +=
					t0_param[5][0] + t0_param[5][1] * event.ssd_channel[2];
			}

			// position
			for (int j = 0; j < 3; ++j) {
				event.be_x[j] = t0.x[be10_indexes[i]][j];
				event.be_y[j] = t0.y[be10_indexes[i]][j];
				event.he_x[j] = t0.x[he4_indexes[i]][j];
				event.he_y[j] = t0.y[he4_indexes[i]][j];
			}

			// fill T0D1 10Be result
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

			// fill T0D1 4He result
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


			// TAF
			// TAF possible 2H number
			int possible_2H = 0;
			if (taf_indexes[i] == -1) {
				event.taf_flag = 1;
				// loop to search for possible 2H
				for (int j = 0; j < 6; ++j) {
					if (taf[j].num == 1 && taf[j].flag[0] == 0x1) {
						++possible_2H;
						taf_indexes[i] = j;
					}
				}
				if (possible_2H > 0 && possible_2H < 4) {
					++possible_2H_num[possible_2H-1];
				} else if (possible_2H >= 4) {
					++possible_2H_num[3];
				}
				++total_possible_2H;
				if (possible_2H != 1) continue;
				// TAFD energy
				event.tafd_energy = taf[taf_indexes[i]].energy[0][0];
				// try to calculate CsI energy
				event.csi_energy = calculator.Energy(
					taf[taf_indexes[i]].theta[0], event.tafd_energy, 150.0
				);
				if (event.csi_energy < 0) continue;
				if (event.csi_energy < 6.0) {
					event.taf_energy = event.tafd_energy + event.csi_energy;
					event.taf_flag = 2;
				} else {
					event.taf_energy = taf[taf_indexes[i]].energy[0][0];
				}
			} else {
				// TAF flag
				event.taf_flag = 0;
				// CsI index
				if (taf[taf_indexes[i]].flag[0] == 0x3) {
					event.csi_index = taf_indexes[i]*2;
				} else if (taf[taf_indexes[i]].flag[0] == 0x5) {
					event.csi_index = taf_indexes[i]*2 + 1;
				} else {
					event.csi_index = -1;
				}
				// CsI channel
				event.csi_channel = taf[taf_indexes[i]].energy[0][1];
				// TAFD energy
				event.tafd_energy = taf[taf_indexes[i]].energy[0][0];
				// TAF-CsI energy
				// calibrated parameters
				double a0 = csi_param[event.csi_index][0];
				double a1 = csi_param[event.csi_index][1];
				double a2 = csi_param[event.csi_index][2];
				// CsI energy
				event.csi_energy = pow(
					(event.csi_channel - a2) / a0,
					1.0 / a1
				);
				// recoil particle energy
				event.taf_energy = event.tafd_energy + event.csi_energy;
			}

			// position
			ROOT::Math::Polar3DVector recoil_position(
				taf[taf_indexes[i]].radius[0],
				taf[taf_indexes[i]].theta[0],
				taf[taf_indexes[i]].phi[0]
			);
			event.d_x = recoil_position.X();
			event.d_y = recoil_position.Y();

			// channel
			event.d_x_channel = tafd[taf_indexes[i]].front_energy[0];
			event.d_y_channel = tafd[taf_indexes[i]].back_energy[0];
			// time
			event.d_x_time = tafd[taf_indexes[i]].front_time[0];
			event.d_y_time = tafd[taf_indexes[i]].back_time[0];
			// strip
			event.d_x_strip = tafd[taf_indexes[i]].front_strip[0];
			event.d_y_strip = tafd[taf_indexes[i]].back_strip[0];

			// PPAC
			// vppac_num = vppac.num;
			xppac_num = xppac.num;
			event.ppac_flag = ppac_flags[i];
			for (int j = 0; j < 3; ++j) {
				event.ppac_x[j] = xppac.x[j] - ppac_correct[0][j];
				event.ppac_y[j] = xppac.y[j] - ppac_correct[1][j];
			}
			// } else if (ppac_flags[i] == 1) {
				// event.tx = vppac.x[3];
			// 	event.ty = vppac.y[3];
			// 	for (int j = 0; j < 3; ++j) {
			// 		event.ppac_x[j] = vppac.x[j];
			// 		event.ppac_y[j] = vppac.y[j];
			// 	}
			// 	event.ppac_xflag = vppac_xflag;
			// 	event.ppac_yflag = vppac_yflag;
			// }
			// PPAC x track
			event.ppac_track[0] = TrackPpac(
				event.ppac_xflag, ppac_xz, event.ppac_x,
				event.t0_energy[0], event.t0_energy[1], event.taf_energy,
				event.be_x[0], event.he_x[0], event.d_x, event.d_y,
				event.tx
			);
			// PPAC y track
			event.ppac_track[1] = TrackPpac(
				event.ppac_yflag, ppac_yz, event.ppac_y,
				event.t0_energy[0], event.t0_energy[1], event.taf_energy,
				event.be_y[0], event.he_y[0], event.d_y, event.d_x,
				event.ty
			);

			// Q
			event.q = ThreeBodyProcess(event);

			// hole
			event.hole[0] = t0.hole[be10_indexes[i]];
			event.hole[1] = t0.hole[he4_indexes[i]];

			// extra information
			// run number
			event.run = run;
			// entry
			event.entry = valid_entries[i];

			opt.Fill();
		}
		t0_tele_file.Close();
	}

	// show possible 2H statistics
	std::cout << "Possible 2H events total " << total_possible_2H << "\n";
	for (int i = 0; i < 4; ++i) {
		std::cout << "With TAF " << i+1 << ": " << possible_2H_num[i] << "\n"; 
	}

	opf.cd();
	// save trees
	opt.Write();
	// close files
	opf.Close();
	return 0;
}