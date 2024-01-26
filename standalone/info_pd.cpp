#include <iostream>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TString.h>
#include <TTree.h>
#include <Math/Vector3D.h>
#include <TGraph.h>
#include <TMultiGraph.h>

#include "include/event/pd_info_event.h"
#include "include/event/channel_event.h"
#include "include/event/ta_event.h"
#include "include/event/t0_event.h"
#include "include/event/particle_event.h"
#include "include/event/dssd_event.h"
#include "include/ppac_track.h"

using namespace ribll;

double PdProcess(
	const PdInfoEvent &event,
	double &c15_kinetic
) {
	// C14 momentum
	double c14_momentum = MomentumFromKinetic(mass_14c, event.t0_energy);
	// C14 momentum vector
	ROOT::Math::XYZVector p_c14(
		event.c_x[0] - event.tx,
		event.c_y[0] - event.ty,
		100.0
	);
	p_c14 = p_c14.Unit() * c14_momentum;

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
	ROOT::Math::XYZVector p_c15 = p_c14 + p_d;

	// 15C momentum
	double c15_momentum = p_c15.R();
	// 15C kinetic energy
	c15_kinetic =
		sqrt(pow(c15_momentum, 2.0) + pow(mass_15c, 2.0)) - mass_15c;

	// Q value
	double q = event.t0_energy + event.taf_energy - c15_kinetic;

	return q;
}


/// @brief search for index from binary flag
/// @param[in] side 0-front, 1-back
/// @param[in] flag binary flag to search
/// @returns valid index of front or back side if only one valid bit is found,
/// 	otherwise, returns -1
int SearchIndex(int side, unsigned short flag) {
	int counts = 0;
	int index = -1;
	int base = side == 0 ? 0x1 : 0x100;
	for (int i = 0; i < 8; ++i) {
		if ((flag & (base << i)) != 0) {
			++counts;
			index = i;
		}
	}
	if (counts != 1) return -1;
	return index;
}


/// @brief get energy from theta
/// @param[in] beam_mass mass of beam
/// @param[in] fragment1_mass mass of fragment1
/// @param[in] fragment2_mass mass of fragment2
/// @param[in] theta angle of fragment1
/// @param[in] q Q value
/// @param[in] beam_kinematic kinematic energy of beam
/// @param[out] fragment1_kinematic1 calculated kinematic energy of fragment1
/// @param[out] fragment1_kinematic2 calculated kinematic energy of fragment2
/// @returns number of roots
///
int EnergyTheta(
	double beam_mass,
	double fragment1_mass,
	double fragment2_mass,
	double theta,
	double q,
	double beam_kinematic,
	double &fragment1_kinematic1,
	double &fragment1_kinematic2
) {
	double t = fragment1_mass+fragment2_mass;
	double r = sqrt(beam_mass*fragment1_mass*beam_kinematic)*cos(theta) / t;
	double s = (
		fragment2_mass*q
		+ fragment2_mass*beam_kinematic
		- beam_mass*beam_kinematic
	) / t;
	double delta = r*r + s;
	if (delta < 0) return 0;
	if (delta == 0) {
		fragment1_kinematic1 = r;
		return 1;
	}
	fragment1_kinematic1 = pow(r + sqrt(delta), 2.0);
	fragment1_kinematic2 = pow(r - sqrt(delta), 2.0);
	if (r - sqrt(delta) < 0) return 1;
	return 2;
}


///	@brief get cos(theta) from energy
/// @param[in] beam_mass mass of beam
/// @param[in] fragment1_mass mass of fragment1
/// @param[in] fragment2_mass mass of fragment2
/// @param[in] q Q value
/// @param[in] beam_kinematic kinematic energy of beam
/// @param[in] fragment1_kinematic kinematic energy of fragment1
/// @returns angle theta of fragment1 and beam
///
double ThetaEnergy(
	double beam_mass,
	double fragment1_mass,
	double fragment2_mass,
	double q,
	double beam_kinematic,
	double fragment1_kinematic
) {
	double fragment2_kinematic = q + beam_kinematic - fragment1_kinematic;

	double momentum_beam = MomentumFromKinetic(beam_mass, beam_kinematic);
	double momentum_fragment1 =
		MomentumFromKinetic(fragment1_mass, fragment1_kinematic);
	double momentum_fragment2 =
		MomentumFromKinetic(fragment2_mass, fragment2_kinematic);

	double cos_theta =
		(
			pow(momentum_fragment1, 2.0)
			+ pow(momentum_beam, 2.0)
			- pow(momentum_fragment2, 2.0)
		)
		/ (2.0 * momentum_fragment1 * momentum_beam);

	return cos_theta;
}


int main() {
	// output file name
	TString output_file_name = TString::Format(
		"%sinfo/pd.root",
		kGenerateDataPath
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// E-theta graph
	TGraph energy_theta[12];
	// theorical E-theta graph
	TGraph theory_energy_theta[2];
	for (int i = 0; i < 2; ++i) {
		theory_energy_theta[i].SetLineColor(kRed);
	}
	// multigraph
	TMultiGraph mg[12];
	// total
	TGraph total_energy_theta;
	TMultiGraph total_mg;
	// output tree
	TTree opt("tree", "information of 15C pd reaction");
	// output event
	PdInfoEvent pd_event;
	// setup output branches
	pd_event.SetupOutput(&opt);

	for (unsigned int run = 428; run <= 452; ++run) {
		if (run == 445) continue;

		std::cout << "Processing run " << run << "\n";
		// input file name
		TString input_file_name = TString::Format(
			"%s%sC15-14C-2H-%04d.root",
			kGenerateDataPath,
			kChannelDir,
			run
		);
		// input file
		TFile input_file(input_file_name, "read");
		// input tree
		TTree *ipt = (TTree*)input_file.Get("tree");
		// input channel event
		ChannelEvent channel;
		// t0 index
		int t0_index;
		// Q value
		double input_q;
		// coincide
		bool coincide;
		// setup input branches
		channel.SetupInput(ipt);
		ipt->SetBranchAddress("t0_index", &t0_index);
		ipt->SetBranchAddress("q", &input_q);
		ipt->SetBranchAddress("coincide", &coincide);

		std::vector<long long> valid_entries;
		std::vector<int> taf_indexes;
		std::vector<int> t0_indexes;
		// std::vector<double> q_values;

		// loop to record information in channel
		for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
			ipt->GetEntry(entry);
			if (!coincide) continue;
			valid_entries.push_back(channel.entry);
			taf_indexes.push_back(channel.taf_index);
			t0_indexes.push_back(t0_index);
		}
		// close file
		input_file.Close();

		// open telescope file to get more information
		// t0 telescope
		TString t0_file_name = TString::Format(
			"%s%st0-telescope-ta-%04u.root",
			kGenerateDataPath,
			kTelescopeDir,
			run
		);
		// t0 file
		TFile t0_file(t0_file_name, "read");
		// t0 tree
		TTree *t0_tree = (TTree*)t0_file.Get("tree");
		// t0 event
		T0Event t0;
		// setup input branches
		t0.SetupInput(t0_tree);

		// TAF particle events
		ParticleEvent tafp[6];
		// add TAF particles as friends
		for (int i = 0; i < 6; ++i) {
			t0_tree->AddFriend(
				TString::Format("taf%dp=tree", i),
				TString::Format(
					"%s%staf%d-particle-ta-%04u.root",
					kGenerateDataPath,
					kParticleDir,
					i,
					run
				)
			);
			// setup input branches
			tafp[i].SetupInput(t0_tree, TString::Format("taf%dp.", i).Data());
		}

		// TAF events
		TaEvent taf[6];
		// add TAF telescopes as friends
		for (int i = 0; i < 6; ++i) {
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
			// setup input branches
			taf[i].SetupInput(t0_tree, TString::Format("taf%d.", i).Data());
		}

		// TAFD event
		DssdFundamentalEvent tafd[6];
		// add TAFD fundamental as friends
		for (int i = 0; i < 6; ++i) {
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
			// setup input branches
			tafd[i].SetupInput(t0_tree, TString::Format("tafd%d.", i).Data());
		}

		// XPPAC events
		ParticleEvent xppac;
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
		// setup input branches
		xppac.SetupInput(t0_tree, "xppac.");
		t0_tree->SetBranchAddress("xppac.xflag", &pd_event.ppac_xflag);
		t0_tree->SetBranchAddress("xppac.yflag", &pd_event.ppac_yflag);


		// T0Dx normalize result friend
		DssdFundamentalEvent norm_events[2];
		// add normalize result friend
		for (int i = 0; i < 2; ++i) {
			t0_tree->AddFriend(
				TString::Format("d%d=tree", i+1),
				TString::Format(
					"%s%st0d%d-result-ta-%04u.root",
					kGenerateDataPath,
					kNormalizeDir,
					i+1,
					run
				)
			);
			norm_events[i].SetupInput(
				t0_tree, TString::Format("d%d.", i+1).Data()
			);
		}

		for (size_t i = 0; i < valid_entries.size(); ++i) {
			t0_tree->GetEntry(valid_entries[i]);

			// TAFD energy
			pd_event.tafd_energy = taf[taf_indexes[i]].energy[0][0];
			// CsI index
			if (taf[taf_indexes[i]].flag[0] == 0x3) {
				pd_event.csi_index = taf_indexes[i]*2;
			} else if (taf[taf_indexes[i]].flag[0] == 0x5) {
				pd_event.csi_index = taf_indexes[i]*2 + 1;
			} else {
				pd_event.csi_index = -1;
			}
			// CsI channel
			pd_event.csi_channel = taf[taf_indexes[i]].energy[0][1];

			// 14C channel
			pd_event.c_channel[0] = t0.energy[t0_indexes[i]][0];
			pd_event.c_channel[1] = t0.energy[t0_indexes[i]][1];

			// calculate energy
			// 14C kinetic energy
			// T0D1 energy
			pd_event.t0_energy =
				t0_param[0][0] + t0_param[0][1] * pd_event.c_channel[0];
			// T0D2 energy
			pd_event.t0_energy +=
				t0_param[1][0] + t0_param[1][1] * pd_event.c_channel[1];

			// CsI calibration parameters
			double a0 = csi_param[pd_event.csi_index][0];
			double a1 = csi_param[pd_event.csi_index][1];
			double a2 = csi_param[pd_event.csi_index][2];
			// 2H kinetic energy
			pd_event.taf_energy = pd_event.tafd_energy + pow(
				(pd_event.csi_channel - a2) / a0,
				1.0 / a1
			);

			// int d1_front_index =
			// 	SearchIndex(0, t0.dssd_flag[t0_indexes[i]][0]);
			// int d1_back_index =
			// 	SearchIndex(1, t0.dssd_flag[t0_indexes[i]][0]);
			// int d2_front_index =
			// 	SearchIndex(0, t0.dssd_flag[t0_indexes[i]][1]);
			// int d2_back_index =
			// 	SearchIndex(1, t0.dssd_flag[t0_indexes[i]][1]);

			// d1x_time = norm_events[0].back_time[d1_back_index];
			// d1y_time = norm_events[0].front_time[d1_front_index];
			// d2x_time = norm_events[1].front_time[d2_front_index];
			// d2y_time = norm_events[1].back_time[d2_back_index];

			// tafd_pi_time = -1e5;
			// for (int j = 0; j < tafd[taf_indexes[i]].front_hit; ++j) {
			// 	if (
			// 		tafd[taf_indexes[i]].front_strip[j]
			// 		== taf[taf_indexes[i]].front_strip[0]
			// 	) {
			// 		tafd_pi_time = tafd[taf_indexes[i]].front_time[j];
			// 	}
			// }
			// tafd_r_time = -1e5;
			// for (int j = 0; j < tafd[taf_indexes[i]].back_hit; ++j) {
			// 	if (
			// 		tafd[taf_indexes[i]].back_strip[j]
			// 		== taf[taf_indexes[i]].back_strip[0]
			// 	) {
			// 		tafd_r_time = tafd[taf_indexes[i]].back_time[j];
			// 	}
			// }

			pd_event.c_x[0] = t0.x[t0_indexes[i]][0];
			pd_event.c_y[0] = t0.y[t0_indexes[i]][0];
			pd_event.c_x[1] = t0.x[t0_indexes[i]][1];
			pd_event.c_y[1] = t0.y[t0_indexes[i]][1];

			pd_event.d_x = tafp[taf_indexes[i]].x[0];
			pd_event.d_y = tafp[taf_indexes[i]].y[0];

			// beam_px = xppac.px[3];
			// beam_py = xppac.py[3];
			// beam_pz = xppac.pz[3];

			// tafd_pi_strip = taf[taf_indexes[i]].front_strip[0];
			// tafd_r_strip = taf[taf_indexes[i]].back_strip[0];

			for (int j = 0; j < 3; ++j) {
				pd_event.ppac_x[j] = xppac.x[j] - ppac_correct[0][j];
				pd_event.ppac_y[j] = xppac.y[j] - ppac_correct[1][j];
			}
			// track PPAC x direction
			pd_event.ppac_x_track = TrackPpac(
				pd_event.ppac_xflag, ppac_xz, pd_event.ppac_x,
				pd_event.t0_energy, pd_event.taf_energy,
				pd_event.c_x[0], pd_event.d_x, pd_event.d_y,
				pd_event.tx
			);
			// track PPAC y direction
			pd_event.ppac_y_track = TrackPpac(
				pd_event.ppac_yflag, ppac_yz, pd_event.ppac_y,
				pd_event.t0_energy, pd_event.taf_energy,
				pd_event.c_y[0], pd_event.d_y, pd_event.d_x,
				pd_event.ty
			);

			pd_event.q = PdProcess(pd_event, pd_event.c15_kinetic);
			pd_event.c15_ex = 1.0064550 - pd_event.q;

			// theta = acos(p1.Dot(p2));
			// // fill energy-theta graph
			// energy_theta[csi_index].AddPoint(theta, taf_energy);
			// total_energy_theta.AddPoint(theta, taf_energy);

			opt.Fill();
		}
		t0_file.Close();
	}


	// // calculate theorical graph
	// for (int i = 0; i < 2; ++i) {
	// 	double q_value = i == 0 ? 1.0064550 : 1.0064550 - 6.0938;
	// 	for (double theta = 0.4; theta < 0.9; theta += 0.01) {
	// 		double kinematic_energy[2];
	// 		int root_num = EnergyTheta(
	// 			mass_15c, mass_2h, mass_14c,
	// 			theta, q_value,
	// 			430.0,
	// 			kinematic_energy[0], kinematic_energy[1]
	// 		);
	// 		if (root_num > 0 && kinematic_energy[0] < 70.0) {
	// 			theory_energy_theta[i].AddPoint(theta, kinematic_energy[0]);
	// 		}
	// 	}
	// 	for (double theta = 0.9; theta > 0.4; theta -= 0.01) {
	// 		double kinematic_energy[2];
	// 		int root_num = EnergyTheta(
	// 			mass_15c, mass_2h, mass_14c,
	// 			theta, q_value,
	// 			430.0,
	// 			kinematic_energy[0], kinematic_energy[1]
	// 		);
	// 		if (root_num == 1) theory_energy_theta[i].AddPoint(theta, kinematic_energy[0]);
	// 		else if (root_num == 2) theory_energy_theta[i].AddPoint(theta, kinematic_energy[1]);
	// 	}
	// }

	for (int i = 0; i < 2; ++i) {
		double q_value = i == 0 ? 1.0064550 : 1.0064550 - 6.0938;
		for (double energy = 0.0; energy < 70.0; energy += 0.5) {
			double cos_theta = ThetaEnergy(
				mass_15c, mass_2h, mass_14c,
				q_value, 430.0, energy
			);
			if (cos_theta > -1.0 && cos_theta < 1.0) {
				theory_energy_theta[i].AddPoint(acos(cos_theta), energy);
			}
		}
	}

	for (int i = 0; i < 12; ++i) {
		mg[i].Add(energy_theta+i, "AP*");
		mg[i].Add(theory_energy_theta, "C");
		mg[i].Add(theory_energy_theta+1, "C");
	}
	total_mg.Add(&total_energy_theta, "AP*");
	total_mg.Add(theory_energy_theta, "C");
	total_mg.Add(theory_energy_theta+1, "C");

	opf.cd();
	theory_energy_theta[0].Write("gt0");
	theory_energy_theta[1].Write("gt1");
	// save energy-theta graph
	for (int i = 0; i < 12; ++i) {
		energy_theta[i].Write(TString::Format("g%d", i));
	}
	for (int i = 0; i < 12; ++i) {
		mg[i].Write(TString::Format("mg%d", i));
	}
	total_energy_theta.Write("gt");
	total_mg.Write("mgt");
	// save trees
	opt.Write();
	// close files
	opf.Close();
	return 0;
}