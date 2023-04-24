#include "include/telescope/t0.h"

#include <TChain.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>
#include <Math/Vector3D.h>

#include "include/event/dssd_event.h"
#include "include/event/ssd_event.h"
#include "include/event/t0_event.h"
#include "include/event/particle_event.h"
#include "include/event/particle_type_event.h"
#include "include/statistics/track_statistics.h"
#include "include/statistics/track_dssd_statistics.h"

namespace ribll {

const double initial_calibration_parameters[12] = {
	0.16, 0.005,
	0.24, 0.0065,
	-0.8, 0.005,
	1.5, 0.0023,
	-0.8, 0.002,
	0.15, 0.0024
};

T0::T0(unsigned int run, const std::string &tag)
: Telescope(run, "t0", tag) {
}


/// @brief track dssd events in T0 telescope
/// @param[in] d1 t0d1 event
/// @param[in] d2 t0d2 event
/// @param[in] d3 t0d3 event
/// @param[in] angle_tolerance angle tolerance
/// @param[out] t0 t0 telescope event
/// @param[out] angle_window array of angle tolerance window
///
void TrackDssdEvent(
	const DssdMergeEvent &d1,
	const DssdMergeEvent &d2,
	const DssdMergeEvent &d3,
	double angle_tolerance,
	T0Event &t0,
	TH1F *angle_window
) {
	// initialize t0 event
	t0.num = 0;
	// fill all d1 events
	for (unsigned short i = 0; i < d1.hit; ++i) {
		t0.layer[t0.num] = 1;
		t0.flag[t0.num] = 0x1;
		t0.energy[t0.num][0] = d1.energy[i];
		t0.x[t0.num][0] = d1.x[i];
		t0.y[t0.num][0] = d1.y[i];
		t0.z[t0.num][0] = d1.z[i];
		++t0.num;
	}
	// fill d2 event in two cases:
	// 1. fill with d1 event if under the angle tolerance
	// 2. fill to empty slot if none match d1 event found
	for (unsigned short i = 0; i < d2.hit; ++i) {
		for (unsigned short j = 0; j < t0.num; ++j) {
			// jump if a d2 event has filled
			if ((t0.flag[j] & 0x2) != 0) continue;
			if ((t0.flag[j] & 0x1) != 0) {
				// found a d1 event, check the angle
				// d1 particle position in lab coordinate
				ROOT::Math::XYZVector d1_pos(
					t0.x[j][0], t0.y[j][0], t0.z[j][0]
				);
				// d2 particle position in lab coordinate
				ROOT::Math::XYZVector d2_pos(d2.x[i], d2.y[i], d2.z[i]);
				double cos_theta =
					d1_pos.Dot(d2_pos) / (d1_pos.R() * d2_pos.R());
				// fill to d1d2 total window
				angle_window[0].Fill(cos_theta);
				if (cos_theta > angle_tolerance) {
					// the angle is acceptable, fill into t0 event
					t0.layer[j]++;
					t0.flag[j] |= 0x2;
					t0.energy[j][1] = d2.energy[i];
					t0.x[j][1] = d2.x[i];
					t0.y[j][1] = d2.y[i];
					t0.z[j][1] = d2.z[i];
					// this d2 event has matched, goto next one
					break;
				}
			} else {
				// none of the d1 events match this d2 event,
				// put the d2 event to empty slot
				t0.layer[j] = 1;
				t0.flag[j] = 0x2;
				t0.energy[j][0] = 0.0;
				t0.energy[j][1] = d2.energy[i];
				t0.x[j][1] = d2.x[i];
				t0.y[j][1] = d2.y[i];
				t0.z[j][1] = d2.z[i];
				++t0.num;
			}
		}
	}
	// fill d3 event in three cases:
	// 1. fill with d1 and d2 events under angle tolerance
	// 2. fill with only d1 event under angle tolerance
	// 3. fill with only d2 event under angle tolerance
	for (unsigned short i = 0; i < d3.hit; ++i) {
		// d3 particle position
		ROOT::Math::XYZVector d3_pos(d3.x[i], d3.y[i], d3.z[i]);
		for (unsigned short j = 0; j < t0.num; ++j) {
			// d1 particle position
			ROOT::Math::XYZVector d1_pos(
				t0.x[j][0], t0.y[j][0], t0.z[j][0]
			);
			// d2 particle position
			ROOT::Math::XYZVector d2_pos(
				t0.x[j][1], t0.y[j][1], t0.z[j][1]
			);
			double cos_theta13 = d1_pos.Dot(d3_pos) / (d1_pos.R()*d3_pos.R());
			double cos_theta23 = d2_pos.Dot(d3_pos) / (d2_pos.R()*d3_pos.R());
			// fill d1d3 angle to total window
			angle_window[1].Fill(cos_theta13);
			// fill d2d3 angle to total window
			angle_window[2].Fill(cos_theta23);
			if (t0.flag[j] == 0x3) {
				// fill d1d3 angle to flag7 window
				angle_window[4].Fill(cos_theta13);
				// fill d2d3 angle to flag7 window
				angle_window[6].Fill(cos_theta23);
				// this slot has d1 and d2 events
				if (
					cos_theta13 > angle_tolerance
					&& cos_theta23 > angle_tolerance
				) {
					// angle check pass, fill d3 event
					t0.layer[j]++;
					t0.flag[j] |= 0x4;
					t0.energy[j][2] = d3.energy[i];
					t0.x[j][2] = d3.x[i];
					t0.y[j][2] = d3.y[i];
					t0.z[j][2] = d3.z[i];
					break;
				}
			} else if (t0.flag[j] == 0x1) {
				// fill d1d3 angle to flag5 window
				angle_window[3].Fill(cos_theta13);
				// this slot only contains d1 event
				if (cos_theta13 > angle_tolerance) {
					// angle check pass, fill d3 event
					t0.layer[j]++;
					t0.flag[j] |= 0x4;
					t0.energy[j][2] = d3.energy[i];
					t0.x[j][2] = d3.x[i];
					t0.y[j][2] = d3.y[i];
					t0.z[j][2] = d3.z[i];
					break;
				}
			} else if (t0.flag[j] == 0x2) {
				// fill d2d3 angle to flag6 window
				angle_window[5].Fill(cos_theta23);
				// this slot only has d2 event
				if (cos_theta23 > angle_tolerance) {
					// angle check pass, fill d3 event
					t0.layer[j]++;
					t0.flag[j] |= 0x4;
					t0.energy[j][2] = d3.energy[i];
					t0.x[j][2] = d3.x[i];
					t0.y[j][2] = d3.y[i];
					t0.z[j][2] = d3.z[i];
					break;
				}
			}
		}
	}
	// discard particles only contains the second layer
	for (unsigned short i = 0; i < t0.num; ++i) {
		if (t0.flag[i] == 0x2) {
			// this particle only contains d2 event,
			// discard it and move the following layer ahead
			for (unsigned short j = i; j < t0.num-1; ++j) {
				t0.layer[j] = t0.layer[j+1];
				t0.flag[j] = t0.flag[j+1];
				for (unsigned short k = 0; k < 3; ++k) {
					t0.energy[j][k] = t0.energy[j+1][k];
					t0.x[j][k] = t0.x[j+1][k];
					t0.y[j][k] = t0.y[j+1][k];
					t0.z[j][k] = t0.z[j+1][k];
				}
				t0.num--;
			}
		}
	}
}


/// @brief track T0 SSD events
/// @param[in] s1 T0S1 event
/// @param[in] s2 T0S2 event
/// @param[in] s3 T0S3 event
/// @param[out] t0 T0 telescope event
void TrackSsdEvent(
	const SsdEvent &s1,
	const SsdEvent &s2,
	const SsdEvent &s3,
	T0Event &t0
) {
	// jump if DSSD not track
	if (t0.num == 0) return;
	t0.ssd_flag = 0;
	if (s1.time > -9e4) {
		t0.ssd_flag |= 0x1;
		t0.ssd_energy[0] = s1.energy;
	}
	if (s2.time > -9e4) {
		t0.ssd_flag |= 0x2;
		t0.ssd_energy[1] = s2.energy;
	}
	if (s3.time > -9e4) {
		t0.ssd_flag |= 0x4;
		t0.ssd_energy[2] = s3.energy;
	}
	return;
}


int T0::Track(double angle_tolerance) {
	// T0D1 merge file name
	TString d1_file_name;
	d1_file_name.Form(
		"%s%st0d1-merge-%s%04u.root",
		kGenerateDataPath,
		kMergeDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// T0D1 file
	TFile d1_file(d1_file_name, "read");
	// tree
	TTree *ipt = (TTree*)d1_file.Get("tree");
	if (!ipt) {
		std::cout << "Error: Get tree from "
			<< d1_file_name << " failed.\n";
		return -1;
	}
	// T0D2 merge file name
	TString d2_file_name;
	d2_file_name.Form(
		"%s%st0d2-merge-%s%04u.root",
		kGenerateDataPath,
		kMergeDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	ipt->AddFriend("d2=tree", d2_file_name);
	// T0D3 merge file name
	TString d3_file_name;
	d3_file_name.Form(
		"%s%st0d3-merge-%s%04u.root",
		kGenerateDataPath,
		kMergeDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	ipt->AddFriend("d3=tree", d3_file_name);
	// T0S1 fundamental file name
	TString s1_file_name;
	s1_file_name.Form(
		"%s%st0s1-fundamental-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	ipt->AddFriend("s1=tree", s1_file_name);
	// T0S2 fundamental file name
	TString s2_file_name;
	s2_file_name.Form(
		"%s%st0s2-fundamental-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	ipt->AddFriend("s2=tree", s2_file_name);
	// T0S3 fundamental file name
	TString s3_file_name;
	s3_file_name.Form(
		"%s%st0s3-fundamental-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	ipt->AddFriend("s3=tree", s3_file_name);
	// d1 event
	DssdMergeEvent d1_event;
	d1_event.SetupInput(ipt);
	// d2 event
	DssdMergeEvent d2_event;
	d2_event.SetupInput(ipt, "d2.");
	// d3 event
	DssdMergeEvent d3_event;
	d3_event.SetupInput(ipt, "d3.");
	// s1 event
	SsdEvent s1_event;
	s1_event.SetupInput(ipt, "s1.");
	// s2 event
	SsdEvent s2_event;
	s2_event.SetupInput(ipt, "s2.");
	// s3 event
	SsdEvent s3_event;
	s3_event.SetupInput(ipt, "s3.");

	// T0 telescope file name
	TString t0_file_name;
	t0_file_name.Form(
		"%s%st0-telescope-%s%04u.root",
		kGenerateDataPath,
		kTelescopeDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// T0 teleescope file
	TFile t0_file(t0_file_name, "recreate");
	// T0 teleescope output tree
	TTree opt("tree", "t0 telescope tree");
	// t0 event
	T0Event t0_event;
	// setup output branches
	t0_event.SetupOutput(&opt);

	// histogram of angle window between DSSD
	TH1F angle_window[]{
		TH1F("ha12", "total angle window of d1d2", 400, 0.6, 1),
		TH1F("ha13", "total angle window of d1d3", 400, 0.6, 1),
		TH1F("ha23", "total angle window of d2d3", 400, 0.6, 1),
		TH1F("ha13a", "angle window of d1d3 under flag 5", 400, 0.6, 1),
		TH1F("ha13b", "angle window of d1d3 under flag 7", 400, 0.6, 1),
		TH1F("ha23a", "angle window of d2d3 under flag 6", 400, 0.6, 1),
		TH1F("ha23b", "angle window of d2d3 under flag 7", 400, 0.6, 1)
	};


	// telescope tracking statistics
	TrackStatistics track_statistics(run_, name_, tag_);
	// DSSD tracking statistics
	TrackDssdStatistics dssd_statistics(run_, name_, tag_);

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Tracking t0   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);
		TrackDssdEvent(
			d1_event, d2_event, d3_event,
			angle_tolerance,
			t0_event,
			angle_window
		);
		TrackSsdEvent(s1_event, s2_event, s3_event, t0_event);
		opt.Fill();
		// update statistics
		if (t0_event.num > 0) {
			// update tracking statistics
			++track_statistics.total;
			if (t0_event.num == 1) ++track_statistics.particle1;
			else if (t0_event.num == 2) ++track_statistics.particle2;
			else if (t0_event.num == 3) ++track_statistics.particle3;
			// update DSSD tracking statistics
			dssd_statistics.total += t0_event.num;
			for (unsigned short i = 0; i < t0_event.num; ++i) {
				if (t0_event.flag[i] == 1) ++dssd_statistics.flag1;
				else if (t0_event.flag[i] == 3) ++dssd_statistics.flag3;
				else if (t0_event.flag[i] == 5) ++dssd_statistics.flag5;
				else if (t0_event.flag[i] == 6) ++dssd_statistics.flag6;
				else if (t0_event.flag[i] == 7) ++dssd_statistics.flag7;
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save tree
	opt.Write();
	// save histograms
	for (auto &hist : angle_window) {
		hist.Write();
	}
	// close files
	t0_file.Close();
	d1_file.Close();

	// save and print statistics
	track_statistics.Write();
	track_statistics.Print();
	dssd_statistics.Write();
	dssd_statistics.Print();

	return 0;
}


int DssdParticleIdentify(
	const T0Event &t0,
	size_t index,
	const std::vector<ParticleCut> &d1d2_cuts,
	const std::vector<ParticleCut> &d2d3_cuts,
	ParticleTypeEvent &type
) {
	if (t0.flag[index] == 0x3) {
		// particle pass 2 DSSD
		for (const auto &cut : d1d2_cuts) {
			if (cut.cut->IsInside(t0.energy[index][1], t0.energy[index][0])) {
				// particle is in cuts, recored the charge and mass number
				type.charge[index] = cut.charge;
				type.mass[index] = cut.mass;
				type.layer[index] = 1;
				return 1;
			}
		}
	} else if (t0.flag[index] == 0x7) {
		// particle pass all 3 DSSD
		for (const auto &cut : d2d3_cuts) {
			if (cut.cut->IsInside(t0.energy[index][2], t0.energy[index][1])) {
				// particle is in cuts, recored the charge and mass number
				type.charge[index] = cut.charge;
				type.mass[index] = cut.mass;
				type.layer[index] = 2;
				return 2;
			}
		}
	}
	type.layer[index] = -1;
	return -1;
}


int SsdParticleIdentify(
	const T0Event &t0,
	size_t index,
	const std::vector<ParticleCut> &d3s1_cuts,
	const std::vector<ParticleCut> &s1s2_cuts,
	const std::vector<ParticleCut> &s2s3_cuts,
	const std::vector<ParticleCut> &s2s3_pass_cuts,
	const std::unique_ptr<TCutG> &s1s2_cross_cut,
	ParticleTypeEvent &type
) {
	// if (t0.flag[index] != 0x7) return -1;
	// the true ssd flag, since there is cross signal in T0S2 from T0S1,
	// check this first and correct the ssd flag
	unsigned short ssd_flag = t0.ssd_flag;
	if (
		ssd_flag == 0x3
		&& s1s2_cross_cut->IsInside(t0.ssd_energy[1], t0.ssd_energy[0])
	) {
		// this is cross signal, correct it
		ssd_flag = 0x1;
	}
	// identify particle in different conditions
	if (ssd_flag == 0x1) {
		// particle stop in the first SSD
		for (const auto &cut : d3s1_cuts) {
			if (cut.cut->IsInside(t0.ssd_energy[0], t0.energy[index][2])) {
				// particle is in cuts, record the charge and mass
				type.charge[index] = cut.charge;
				type.mass[index] = cut.mass;
				type.layer[index] = 3;
				return 3;
			}
		}
	} else if (ssd_flag == 0x3) {
		// particle stop in the second SSD
		for (const auto &cut : s1s2_cuts) {
			if (cut.cut->IsInside(t0.ssd_energy[1], t0.ssd_energy[0])) {
				// particle is in cuts, record the charge and mass
				type.charge[index] = cut.charge;
				type.mass[index] = cut.mass;
				type.layer[index] = 4;
				return 4;
			}
		}
	} else if (ssd_flag == 0x7) {
		// particle stop in the third SSD or CsI(Tl)
		for (const auto &cut : s2s3_cuts) {
			if (cut.cut->IsInside(t0.ssd_energy[2], t0.ssd_energy[1])) {
				type.charge[index] = cut.charge;
				type.mass[index] = cut.mass;
				type.layer[index] = 5;
				return 5;
			}
		}
		for (const auto &cut : s2s3_pass_cuts) {
			if (cut.cut->IsInside(t0.ssd_energy[2], t0.ssd_energy[1])) {
				type.charge[index] = cut.charge;
				type.mass[index] = cut.mass;
				type.layer[index] = 6;
				return 6;
			}
		}
	}
	type.layer[index] = -1;
	return -1;
}


int T0::ParticleIdentify() {
	// input t0 telescope file name
	TString t0_file_name;
	t0_file_name.Form(
		"%s%st0-telescope-%s%04u.root",
		kGenerateDataPath,
		kTelescopeDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// t0 file
	TFile t0_file(t0_file_name, "read");
	// t0 tree
	TTree *ipt = (TTree*)t0_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< t0_file_name << " failed.\n";
		t0_file.Close();
		return -1;
	}
	// input event
	T0Event t0_event;
	// setup input branches
	t0_event.SetupInput(ipt);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%st0-particle-type-%s%04u.root",
		kGenerateDataPath,
		kParticleIdentifyDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "particle type");
	// output data
	ParticleTypeEvent type_event;
	// setup output branches
	type_event.SetupOutput(&opt);

	// T0D1-D2 cuts
	std::vector<ParticleCut> d1d2_cuts;
	d1d2_cuts.push_back({2, 4, ReadCut("d1d2", "4He")});
	d1d2_cuts.push_back({3, 6, ReadCut("d1d2", "6Li")});
	d1d2_cuts.push_back({3, 7, ReadCut("d1d2", "7Li")});
	d1d2_cuts.push_back({4, 7, ReadCut("d1d2", "7Be")});
	d1d2_cuts.push_back({4, 9, ReadCut("d1d2", "9Be")});
	d1d2_cuts.push_back({4, 10, ReadCut("d1d2", "10Be")});
	d1d2_cuts.push_back({5, 10, ReadCut("d1d2", "10B")});
	d1d2_cuts.push_back({5, 11, ReadCut("d1d2", "11B")});
	d1d2_cuts.push_back({5, 12, ReadCut("d1d2", "12B")});
	d1d2_cuts.push_back({5, 13, ReadCut("d1d2", "13B")});
	d1d2_cuts.push_back({6, 12, ReadCut("d1d2", "12C")});
	d1d2_cuts.push_back({6, 13, ReadCut("d1d2", "13C")});
	d1d2_cuts.push_back({6, 14, ReadCut("d1d2", "14C")});
	// T0D2-D3 cuts
	std::vector<ParticleCut> d2d3_cuts;
	d2d3_cuts.push_back({2, 4, ReadCut("d2d3", "4He")});
	d2d3_cuts.push_back({3, 7, ReadCut("d2d3", "7Li")});
	d2d3_cuts.push_back({4, 10, ReadCut("d2d3", "10Be")});
	// T0D3-S1 cuts
	std::vector<ParticleCut> d3s1_cuts;
	d3s1_cuts.push_back({2, 4, ReadCut("d3s1", "4He")});

	// T0S1-S2 cuts
	std::vector<ParticleCut> s1s2_cuts;
	s1s2_cuts.push_back({2, 4, ReadCut("s1s2", "4He")});
	// T0S2-S3 cuts
	std::vector<ParticleCut> s2s3_cuts;
	s2s3_cuts.push_back({2, 4, ReadCut("s2s3", "4He")});
	// T0S2-S3 pass cuts
	std::vector<ParticleCut> s2s3_pass_cuts;
	s2s3_pass_cuts.push_back({2, 4, ReadCut("s2s3-t", "4He")});
	// special cut of S1-S2 interaction
	std::unique_ptr<TCutG> s1s2_cross_cut{ReadCut("s1s2", "cross")};

	//statistics
	long long total = 0;
	long long id11 = 0;
	long long id21 = 0;
	long long id22 = 0;

	// total number of particles
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Identifying T0 particles   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);
		// initialize output event
		type_event.num = t0_event.num;
		for (unsigned short i = 0; i < type_event.num; ++i) {
			type_event.charge[i] = 0;
			type_event.mass[i] = 0;
			type_event.layer[i] = -1;
		}
		// identify particles
		if (t0_event.num == 1) {
			++total;
			int identify = DssdParticleIdentify(
				t0_event, 0, d1d2_cuts, d2d3_cuts, type_event
			);
			if (
				identify <= 0
				&& t0_event.flag[0] == 0x7
				&& t0_event.ssd_flag != 0
			) {
				identify = SsdParticleIdentify(
					t0_event, 0,
					d3s1_cuts, s1s2_cuts, s2s3_cuts, s2s3_pass_cuts,
					s1s2_cross_cut,
					type_event
				);
			}
			if (identify > 0) ++id11;
		} else if (t0_event.num == 2) {
			++total;
			int identify0 = DssdParticleIdentify(
				t0_event, 0, d1d2_cuts, d2d3_cuts, type_event
			);
			int identify1 = DssdParticleIdentify(
				t0_event, 1, d1d2_cuts, d2d3_cuts, type_event
			);
			if (identify0 <= 0 && identify1 > 0 && t0_event.flag[0] == 0x7) {
				int identify = SsdParticleIdentify(
					t0_event, 0,
					d3s1_cuts, s1s2_cuts, s2s3_cuts, s2s3_pass_cuts,
					s1s2_cross_cut,
					type_event
				);
				if (identify > 0) ++id22;
				else ++id21;
			} else if (
				identify0 > 0 && identify1 <= 0
				&& t0_event.flag[1] == 0x7
			) {
				int identify = SsdParticleIdentify(
					t0_event, 1,
					d3s1_cuts, s1s2_cuts, s2s3_cuts, s2s3_pass_cuts,
					s1s2_cross_cut,
					type_event
				);
				if (identify > 0) ++id22;
				else ++id21;
			} else if (identify0 > 0 && identify1 > 0) {
				++id22;
			}
		}
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// write tree
	opt.Write();
	// close files
	opf.Close();
	t0_file.Close();

	// show statistics
	std::cout << "Particle identify in T0\n"
		<< "Idnetify 1/1 " << id11 << " / " << total
		<< "  " << double(id11) / double(total) << "\n"
		<< "Idnetify 1/2 " << id21 << " / " << total
		<< "  " << double(id21) / double(total) << "\n"
		<< "Idnetify 2/2 " << id22 << " / " << total
		<< "  " << double(id22) / double(total) << "\n";
	return 0;
}



struct ParticlePidInfo {
	// layer, 0 for d1d2, 1 for d2d3, 2 for d3s1, 3 for s1s2, 4 for s2s3
	unsigned short layer;
	// particle type, 0 for 1H, 1 for 2H, 2 for 4He,
	// 3 for 7Li, 4 for 10Be, 5 for 14C
	unsigned short type;
	// charge number of this particle
	unsigned short charge;
	// mass number of this particle
	unsigned short mass;
	// left bound of pid
	double left;
	// right bound of pid
	double right;
	// offset of this layer
	double offset;
};


const std::vector<ParticlePidInfo> pid_info {
	// d1d2
	{0, 2, 2, 4, 4000.0, 8000.0, 0.0},
	{0, 3, 3, 7, 8000.0, 13'000.0, 0.0},
	{0, 4, 4, 10, 13'000.0, 24'000.0, 0.0},
	{0, 5, 6, 14, 36'000.0, 52'000.0, 0.0},
	// d2d3
	{1, 2, 2, 4, 5000.0, 9000.0, 55'000.0},
	{1, 3, 3, 7, 10'000.0, 16'000.0, 55'000.0},
	{1, 4, 4, 10, 20'000.0, 26'000.0, 55'000.0},
	// d3s1
	{2, 1, 1, 2, 1900.0, 2400.0, 85'000.0},
	{2, 2, 2, 4, 5200.0, 9500.0, 85'000.0},
	// s1s2
	{3, 0, 1, 1, 3500.0, 5000.0, 95'000.0},
	{3, 2, 2, 4, 13'500.0, 25'000.0, 95'000.0},
	// s2s3
	{4, 0, 1, 1, 3500.0, 5500.0, 120'000.0},
	{4, 2, 2, 4, 13'500.0, 22'500.0, 120'000.0}
};


class PidFitFunc {
public:
	PidFitFunc(
		const std::string &telescope,
		const std::vector<std::string> &projectiles
	) {
		for (const std::string &p : projectiles) {
			calculators_.emplace_back(telescope, p);
		}
	}

	double operator()(double *x, double *par) const {
		// identify particle
		for (const auto &info : pid_info) {
			if (
				x[0] > info.left + info.offset
				&& x[0] < info.right + info.offset
			) {
				double de = par[info.layer*2] + par[info.layer*2+1] * (x[0] - info.offset);
				double e = calculators_[info.type].Energy(info.layer, de);
				return (e - par[info.layer*2+2]) / par[info.layer*2+3];
			}
		}
		return 0.0;
	}
private:
	std::vector<elc::DeltaEnergyCalculator> calculators_;
};


int T0::Calibrate(unsigned int end_run) {
	TChain ipt("chain", "chain of T0 events");
	for (unsigned int i = run_; i <= end_run; ++i) {
		if (i == 628) continue;
		ipt.AddFile(TString::Format(
			"%s%st0-telescope-%s%04u.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			tag_.empty() ? "" : (tag_+"-").c_str(),
			i
		));
	}

	// input telescope event
	T0Event t0_event;
	// setup input branches
	t0_event.SetupInput(&ipt);

	// output calibration root file name
	TString calibration_file_name;
	calibration_file_name.Form(
		"%s%st0-calibration-%s%04u-%04u.root",
		kGenerateDataPath,
		kCalibrationDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_,
		end_run
	);
	// output calibration file
	TFile calibration_file(calibration_file_name, "recreate");
	// output E-VS-dE graph, with offset for different layers in dE
	TGraph e_vs_de_offset;

	// T0D1-D2 cuts
	std::vector<ParticleCut> d1d2_cuts;
	d1d2_cuts.push_back({2, 4, ReadCut("d1d2", "4He")});
	d1d2_cuts.push_back({3, 7, ReadCut("d1d2", "7Li")});
	d1d2_cuts.push_back({4, 10, ReadCut("d1d2", "10Be")});
	// d1d2_cuts.push_back({5, 12, ReadCut("d1d2", "12B")});
	d1d2_cuts.push_back({6, 14, ReadCut("d1d2", "14C")});
	// T0D2-D3 cuts
	std::vector<ParticleCut> d2d3_cuts;
	d2d3_cuts.push_back({2, 4, ReadCut("d2d3", "4He")});
	d2d3_cuts.push_back({3, 7, ReadCut("d2d3", "7Li")});
	d2d3_cuts.push_back({4, 10, ReadCut("d2d3", "10Be")});
	// T0D3-S1 cuts
	std::vector<ParticleCut> d3s1_cuts;
	d3s1_cuts.push_back({1, 2, ReadCut("d3s1", "2H")});
	d3s1_cuts.push_back({2, 4, ReadCut("d3s1", "4He")});

	// T0S1-S2 cuts
	std::vector<ParticleCut> s1s2_cuts;
	s1s2_cuts.push_back({1, 1, ReadCut("s1s2", "1H")});
	s1s2_cuts.push_back({2, 4, ReadCut("s1s2", "4He")});
	// T0S2-S3 cuts
	std::vector<ParticleCut> s2s3_cuts;
	s2s3_cuts.push_back({1, 1, ReadCut("s2s3", "1H")});
	s2s3_cuts.push_back({2, 4, ReadCut("s2s3", "4He")});
	// T0S2-S3 pass cuts
	std::vector<ParticleCut> s2s3_pass_cuts;
	// special cut of S1-S2 interaction
	std::unique_ptr<TCutG> s1s2_cross_cut{ReadCut("s1s2", "cross")};

	// type event for filter
	ParticleTypeEvent type_event;
	type_event.num = 1;

	// total number of entries
	long long entries = ipt.GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling events to graph   0%%");
	fflush(stdout);
	// loop to fill events to graph
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt.GetEntry(entry);
		// initialize type event
		type_event.charge[0] = type_event.mass[0] = 0;
		// only use num==1 events to calibrate
		if (t0_event.num != 1) continue;
		int layer = DssdParticleIdentify(
			t0_event, 0, d1d2_cuts, d2d3_cuts, type_event
		);
		double delta_energy = -1e5;
		double energy = -1e5;
		// check layer and get dE and E
		if (layer == 1 || layer == 2) {
			delta_energy = t0_event.energy[0][layer-1];
			energy = t0_event.energy[0][layer];
			for (const auto &info : pid_info) {
				if (
					info.layer == layer - 1
					&& info.charge == type_event.charge[0]
					&& info.mass == type_event.mass[0]
					&& delta_energy > info.left
					&& delta_energy < info.right
				) {
					e_vs_de_offset.AddPoint(
						delta_energy + info.offset, energy
					);
					break;
				}
			}
		}
		layer = SsdParticleIdentify(
			t0_event, 0,
			d3s1_cuts, s1s2_cuts, s2s3_cuts, s2s3_pass_cuts, s1s2_cross_cut,
			type_event
		);
		if (layer == 3) {
			delta_energy = t0_event.energy[0][2];
			energy = t0_event.ssd_energy[0];
		} else if (layer == 4 || layer == 5) {
			delta_energy = t0_event.ssd_energy[layer-4];
			energy = t0_event.ssd_energy[layer-3];
		}
		if (layer < 3 || layer > 5) continue;
		for (const auto &info : pid_info) {
			if (
				info.layer == layer - 1
				&& info.charge == type_event.charge[0]
				&& info.mass == type_event.mass[0]
				&& delta_energy > info.left
				&& delta_energy < info.right
			) {
				e_vs_de_offset.AddPoint(
					delta_energy + info.offset, energy
				);
				break;
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");


	std::vector<std::string> projectiles{
		"1H", "2H", "4He", "7Li", "10Be", "14C"
	};
	PidFitFunc pid_fit("t0", projectiles);
	// fit function
	TF1 fcali("fcali", pid_fit, 0.0, 150'000.0, 12);
	fcali.SetNpx(10000);
	for (size_t i = 0; i < 12; ++i) {
		fcali.SetParameter(i, initial_calibration_parameters[i]);
	}
	// fit
	e_vs_de_offset.Fit(&fcali, "R+");

	// get parameters from fitting
	for (size_t i = 0; i < 12; ++i) {
		cali_params_[i] = fcali.GetParameter(i);
	}
	// save parameters to txt file
	if (WriteCalibrateParameters()) {
		calibration_file.Close();
		// telescope_file.Close();
		return -1;
	}

	// save graph
	calibration_file.cd();
	e_vs_de_offset.Write("gcali");
	// close files
	calibration_file.Close();
	// telescope_file.Close();

	return 0;
}


int T0::Rebuild() {
	// telescope file name
	TString telescope_file_name;
	telescope_file_name.Form(
		"%s%st0-telescope-%s%04u.root",
		kGenerateDataPath,
		kTelescopeDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// telescope file
	TFile telescope_file(telescope_file_name, "read");
	// input telescope tree
	TTree *ipt = (TTree*)telescope_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< telescope_file_name << " failed.\n";
		return -1;
	}
	// particle type file name
	TString particle_type_file_name;
	particle_type_file_name.Form(
		"%s%st0-particle-type-%s%04u.root",
		kGenerateDataPath,
		kParticleIdentifyDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// add friend
	ipt->AddFriend("type=tree", particle_type_file_name);
	// input telescope event
	T0Event t0_event;
	// input type event
	ParticleTypeEvent type_event;
	// setup input branches
	t0_event.SetupInput(ipt);
	type_event.SetupInput(ipt, "type.");

	// output file name
	TString particle_file_name;
	particle_file_name.Form(
		"%s%st0-particle-%s%04u.root",
		kGenerateDataPath,
		kParticleDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile particle_file(particle_file_name, "recreate");
	// particle tree
	TTree opt("tree", "t0 particles");
	// output particle event
	ParticleEvent particle_event;
	// setup output branches
	particle_event.SetupOutput(&opt);

	// read calibrate parameters
	if (ReadCalibrateParameters()) {
		std::cerr << "Errro: Read calibrate parameters failed.\n";
		return -1;
	}

	// CsI energy calculator
	elc::CsiEnergyCalculator csi_calculator("4He");

	// totla number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Rebuilding particle   0%%");
	fflush(stdout);
	// loop events
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);
		// initialize particle event
		particle_event.num = 0;
		for (unsigned short i = 0; i < t0_event.num; ++i) {
			// jump confuesd particles
			if (type_event.charge[i] <= 0 || type_event.mass[i] <= 0) continue;

			// fill charge number
			particle_event.charge[particle_event.num] = type_event.charge[i];
			// fill mass number
			particle_event.mass[particle_event.num] = type_event.mass[i];

			// calculate energy from all Si detectors
			double energy =
				TotalEnergy(t0_event, type_event, i, csi_calculator);
			// fill energy
			particle_event.energy[particle_event.num] = energy;

			// set particle position, T0D1 may has the best resolution
			particle_event.x[particle_event.num] = t0_event.x[i][0];
			particle_event.y[particle_event.num] = t0_event.y[i][0];
			particle_event.z[particle_event.num] = t0_event.z[i][0];

			// if (t0_event.flag[i] == 0x3) {
			// 	particle_event.x[particle_event.num] = t0_event.x[i][1];
			// 	particle_event.y[particle_event.num] = t0_event.y[i][1];
			// 	particle_event.z[particle_event.num] = t0_event.z[i][1];
			// } else if (t0_event.flag[i] == 0x7) {
			// 	particle_event.x[particle_event.num] = t0_event.x[i][2];
			// 	particle_event.y[particle_event.num] = t0_event.y[i][2];
			// 	particle_event.z[particle_event.num] = t0_event.z[i][2];
			// } else {
			// 	std::cerr << "Error: Should not be here in T0::Particle.\n";
			// 	return -1;
			// }


			// leave momentum empty
			particle_event.px[particle_event.num] = 0.0;
			particle_event.py[particle_event.num] = 0.0;
			particle_event.pz[particle_event.num] = 0.0;

			// // calculate momentum value from energy
			// double momentum = MomentumFromEnergy(energy, AccurateMass(type_event.mass[i]));
			// // get momentum direction from Si strips
			// ROOT::Math::XYZVector direction(0.0, 0.0, 0.0);
			// if (t0_event.flag[i] == 0x3) {
			// 	// direction.SetXYZ(
			// 	// 	t0_event.x[i][1] - t0_event.x[i][0],
			// 	// 	t0_event.y[i][1] - t0_event.y[i][0],
			// 	// 	t0_event.z[i][1] - t0_event.z[i][0]
			// 	// );
			// 	direction.SetXYZ(
			// 		t0_event.x[i][1], t0_event.y[i][1], t0_event.z[i][1]
			// 	);
			// } else if (t0_event.flag[i] == 0x7) {
			// 	// direction.SetXYZ(
			// 	// 	t0_event.x[i][2] - t0_event.x[i][0],
			// 	// 	t0_event.y[i][2] - t0_event.y[i][0],
			// 	// 	t0_event.z[i][2] - t0_event.z[i][0]
			// 	// );
			// 	direction.SetXYZ(
			// 		t0_event.x[i][2], t0_event.y[i][2], t0_event.z[i][2]
			// 	);
			// } else {
			// 	std::cerr << "Error: Should not be here in T0::Particle.\n";
			// 	return -1;
			// }
			// direction = direction.Unit();
			// // fill px, py, pz
			// particle_event.px[particle_event.num] = direction.X() * momentum;
			// particle_event.py[particle_event.num] = direction.Y() * momentum;
			// particle_event.pz[particle_event.num] = direction.Z() * momentum;
			++particle_event.num;
		}
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save particle tree
	particle_file.cd();
	opt.Write();
	// close files
	particle_file.Close();
	telescope_file.Close();

	return 0;
}


double T0::TotalEnergy(
	const T0Event &t0,
	const ParticleTypeEvent &type,
	const size_t index,
	const elc::CsiEnergyCalculator &csi_calculator
) const {
	// total energy
	double energy = 0.0;
	// DSSD layers particle hit
	size_t dssd_layer = type.layer[index] == 1 ? 1 : 2;
	for (size_t i = 0; i <= dssd_layer; ++i) {
		// add calibrated DSSD energy
		energy += CaliEnergy(i, t0.energy[index][i]);
	}
	// return energy if particle stop at T0D2
	if (dssd_layer == 1) return energy;

	size_t ssd_layer = type.layer[index] < 5 ? type.layer[index] : 5;
	ssd_layer -= 2;
	for (size_t i = 0; i < ssd_layer; ++i) {
		// add calibrated SSD energy
		energy += CaliEnergy(i+3, t0.ssd_energy[i]);
	}
	// return energy if particle stop in Si
	if (type.layer[index] != 6) return energy;

	double thickness = 0.0;
	for (size_t i = 0; i < 6; ++i) thickness += t0_thickness[i];
	// calculate energy in CsI(Tl)
	double csi_energy = csi_calculator.Energy(
		0.0, energy, thickness
	);

	// add calculated CsI energy
	energy += csi_energy;
	// return energy
	return energy;
}


}		// namespace ribll