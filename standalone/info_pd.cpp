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

#include "include/event/channel_event.h"
#include "include/event/ta_event.h"
#include "include/event/t0_event.h"
#include "include/event/particle_event.h"
#include "include/event/dssd_event.h"

using namespace ribll;


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


/// @brief get momentum from kinematic energy
/// @param[in] mass mass of particle
/// @param[in] kinematic kinematic energy of particle
/// @returns momentum of particle
///
inline double MomentumFromKinematic(double mass, double kinematic) {
	return sqrt((2.0 * mass + kinematic) * kinematic);
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

	double momentum_beam = MomentumFromKinematic(beam_mass, beam_kinematic);
	double momentum_fragment1 =
		MomentumFromKinematic(fragment1_mass, fragment1_kinematic);
	double momentum_fragment2 =
		MomentumFromKinematic(fragment2_mass, fragment2_kinematic);

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
	// output data
	// indexes
	int csi_index;
	int layer;
	// energy
	double taf_energy;
	double tafd_energy;
	// double csi_energy;
	// double t0_energy;
	// double d1_energy;
	// double d2_energy;
	// channel
	double csi_channel;
	double d1_channel;
	double d2_channel;
	// double d1x_channel;
	// double d1y_channel;
	// double d2x_channel;
	// double d2y_channel;
	// time
	double d1x_time;
	double d1y_time;
	double d2x_time;
	double d2y_time;
	double tafd_pi_time;
	double tafd_r_time;
	// position
	double d1x;
	double d1y;
	double d2x;
	double d2y;
	// recoil position
	double rx;
	double ry;
	// target position
	double tx;
	double ty;
	// beam direction
	double beam_px;
	double beam_py;
	double beam_pz;
	// strips
	// unsigned short d1x_strip;
	// unsigned short d1y_strip;
	// unsigned short d2x_strip;
	// unsigned short d2y_strip;
	unsigned short tafd_pi_strip;
	unsigned short tafd_r_strip;
	// bool d1x_merge;
	// bool d1y_merge;
	// bool d2x_merge;
	// bool d2y_merge;
	// bool tafd_r_merge;
	// bool tafd_pi_merge;
	unsigned short ppac_xflag;
	unsigned short ppac_yflag;
	double ppac_x[3];
	double ppac_y[3];
	// state
	int state;
	double q_value;
	// theta
	double theta;

	// setup output branches
	// indexes and layers
	opt.Branch("csi_index", &csi_index, "ci/I");
	opt.Branch("layer", layer, "layer/I");
	// energy
	opt.Branch("taf_energy", &taf_energy, "tafe/D");
	opt.Branch("tafd_energy", &tafd_energy, "tafde/D");
	// opt.Branch("csi_energy", &csi_energy, "csie/D");
	// opt.Branch("t0_energy", &t0_energy, "t0e/D");
	// opt.Branch("d1_energy", &d1_energy, "d1e/D");
	// opt.Branch("d2_energy", &d2_energy, "d2e/D");
	// channel
	opt.Branch("csi_channel", &csi_channel, "csic/D");
	opt.Branch("d1_channel", &d1_channel, "d1c/D");
	opt.Branch("d2_channel", &d2_channel, "d2c/D");
	// opt.Branch("d1_x_channel", &d1x_channel, "d1xc/D");
	// opt.Branch("d1_y_channel", &d1y_channel, "d1yc/D");
	// opt.Branch("d2_x_channel", &d2x_channel, "d2xc/D");
	// opt.Branch("d2_y_channel", &d2y_channel, "d2yc/D");
	// time
	opt.Branch("d1_x_time", &d1x_time, "d1xt/D");
	opt.Branch("d1_y_time", &d1y_time, "d1yt/D");
	opt.Branch("d2_x_time", &d2x_time, "d2xt/D");
	opt.Branch("d2_y_time", &d2y_time, "d2yt/D");
	opt.Branch("tafd_pi_time", &tafd_pi_time, "tafdpt/D");
	opt.Branch("tafd_r_time", &tafd_r_time, "tafdrt/D");
	// position
	opt.Branch("d1_x", &d1x, "d1x/D");
	opt.Branch("d1_y", &d1y, "d1y/D");
	opt.Branch("d2_x", &d2x, "d2x/D");
	opt.Branch("d2_y", &d2y, "d2y/D");
	opt.Branch("rx", &rx, "rx/D");
	opt.Branch("ry", &ry, "ry/D");
	opt.Branch("tx", &tx, "tx/D");
	opt.Branch("ty", &ty, "ty/D");
	opt.Branch("beam_px", &beam_px, "bpx/D");
	opt.Branch("beam_py", &beam_py, "bpy/D");
	opt.Branch("beam_pz", &beam_pz, "bpz/D");
	// strips
	// opt.Branch("d1x_strip", &d1x_strip, "d1xs/s");
	// opt.Branch("d1y_strip", &d1y_strip, "d1ys/s");
	// opt.Branch("d2x_strip", &d2x_strip, "d2xs/s");
	// opt.Branch("d2y_strip", &d2y_strip, "d2ys/s");
	opt.Branch("tafd_pi_strip", &tafd_pi_strip, "tafdps/s");
	opt.Branch("tafd_r_strip", &tafd_r_strip, "tafdrs/s");
	// opt.Branch("d1x_merge", &d1x_merge, "d1xm/O");
	// opt.Branch("d1y_merge", &d1y_merge, "d1ym/O");
	// opt.Branch("d2x_merge", &d2x_merge, "d2xm/O");
	// opt.Branch("d2y_merge", &d2y_merge, "d2ym/O");
	// opt.Branch("tafd_r_merge", &tafd_r_merge, "tafdrm/O");
	// opt.Branch("tafd_pi_merge", &tafd_pi_merge, "tafdpm/O");
	opt.Branch("ppac_xflag", &ppac_xflag, "pxflag/s");
	opt.Branch("ppac_yflag", &ppac_yflag, "pyflag/s");
	opt.Branch("ppac_x", ppac_x, "ppacx[3]/D");
	opt.Branch("ppac_y", ppac_y, "ppacy[3]/D");
	// state
	opt.Branch("state", &state, "state/I");
	opt.Branch("q", &q_value, "q/D");
	// theta
	opt.Branch("theta", &theta, "theta/D");

	for (unsigned int run = 421; run <= 452; ++run) {
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
		std::vector<double> q_values;

		// loop to record information in channel
		for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
			ipt->GetEntry(entry);
			if (!coincide) continue;
			valid_entries.push_back(channel.entry);
			taf_indexes.push_back(channel.taf_index);
			t0_indexes.push_back(t0_index);
			// double q = channel.parent_energy + 1.006
			// 	- channel.daughter_energy[0] - channel.daughter_energy[1];
			q_values.push_back(input_q);
			// if (q > -3.0 && q < 2.5) states.push_back(0);
			// else if (q > 3.0 && q < 9.0) states.push_back(1);
			// else states.push_back(-1);
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
		t0_tree->SetBranchAddress("xppac.xflag", &ppac_xflag);
		t0_tree->SetBranchAddress("xppac.yflag", &ppac_yflag);


		// T0Dx normalize result friend
		DssdFundamentalEvent norm_events[2];
		// add normalize friend
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

			taf_energy = tafp[taf_indexes[i]].energy[0];
			tafd_energy = taf[taf_indexes[i]].energy[0][0];
			if (taf[taf_indexes[i]].flag[0] == 0x3) {
				csi_index = taf_indexes[i]*2;
			} else if (taf[taf_indexes[i]].flag[0] == 0x5) {
				csi_index = taf_indexes[i]*2 + 1;
			} else {
				csi_index = -1;
			}
			layer = 0;

			csi_channel = taf[taf_indexes[i]].energy[0][1];
			d1_channel = t0.energy[t0_indexes[i]][0];
			d2_channel = t0.energy[t0_indexes[i]][1];

			// taf_energy = tafd_energy + pow((csi_channel-200.0)/200.0, 1.0/0.97);

			int d1_front_index =
				SearchIndex(0, t0.dssd_flag[t0_indexes[i]][0]);
			int d1_back_index =
				SearchIndex(1, t0.dssd_flag[t0_indexes[i]][0]);
			int d2_front_index =
				SearchIndex(0, t0.dssd_flag[t0_indexes[i]][1]);
			int d2_back_index =
				SearchIndex(1, t0.dssd_flag[t0_indexes[i]][1]);

			d1x_time = norm_events[0].back_time[d1_back_index];
			d1y_time = norm_events[0].front_time[d1_front_index];
			d2x_time = norm_events[1].front_time[d2_front_index];
			d2y_time = norm_events[1].back_time[d2_back_index];

			tafd_pi_time = -1e5;
			for (int j = 0; j < tafd[taf_indexes[i]].front_hit; ++j) {
				if (
					tafd[taf_indexes[i]].front_strip[j]
					== taf[taf_indexes[i]].front_strip[0]
				) {
					tafd_pi_time = tafd[taf_indexes[i]].front_time[j];
				}
			}
			tafd_r_time = -1e5;
			for (int j = 0; j < tafd[taf_indexes[i]].back_hit; ++j) {
				if (
					tafd[taf_indexes[i]].back_strip[j]
					== taf[taf_indexes[i]].back_strip[0]
				) {
					tafd_r_time = tafd[taf_indexes[i]].back_time[j];
				}
			}

			d1x = t0.x[t0_indexes[i]][0];
			d1y = t0.y[t0_indexes[i]][0];
			d2x = t0.x[t0_indexes[i]][1];
			d2y = t0.y[t0_indexes[i]][1];

			rx = tafp[taf_indexes[i]].x[0];
			ry = tafp[taf_indexes[i]].y[0];

			tx = xppac.x[3];
			ty = xppac.y[3];

			beam_px = xppac.px[3];
			beam_py = xppac.py[3];
			beam_pz = xppac.pz[3];

			tafd_pi_strip = taf[taf_indexes[i]].front_strip[0];
			tafd_r_strip = taf[taf_indexes[i]].back_strip[0];

			for (int i = 0; i < 3; ++i) {
				ppac_x[i] = xppac.x[i];
				ppac_y[i] = xppac.y[i];
			}

			q_value = q_values[i];

			if (q_value < 3.5) state = 0;
			else state = 1;

			ROOT::Math::XYZVector p1(rx-tx, ry-ty, 135.0);
			p1 = p1.Unit();
			ROOT::Math::XYZVector p2(xppac.px[3], xppac.py[3], xppac.pz[3]);
			p2 = p2.Unit();
			theta = acos(p1.Dot(p2));
			// fill energy-theta graph
			energy_theta[csi_index].AddPoint(theta, taf_energy);
			total_energy_theta.AddPoint(theta, taf_energy);

			opt.Fill();
		}
		t0_file.Close();
	}


	constexpr double u = 931.494;
	// constexpr double mass_1h = u * 1.0072764520;
	constexpr double mass_2h = u * 2.0135531980;
	constexpr double mass_14c = u * 13.9999505089;
	constexpr double mass_15c = u * 15.0073077289;
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