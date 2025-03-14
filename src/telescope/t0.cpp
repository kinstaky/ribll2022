#include "include/telescope/t0.h"

#include <queue>

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
	-0.1, 0.007,
	0.5, 0.007,
	-0.1, 0.0055,
	0.05, 0.0023,
	0.0, 0.0025,
	0.0, 0.0025
};

T0::T0(unsigned int run, const std::string &tag)
: Telescope(run, "t0", tag) {
}


/// @brief track dssd events in T0 telescope
/// @param[in] d1 t0d1 event
/// @param[in] d2 t0d2 event
/// @param[in] d3 t0d3 event
/// @param[in] hole hole flag
/// @param[out] t0 t0 telescope event
/// @param[out] offset_window array of offset tolerance window
///
void TrackDssdEvent(
	const DssdMergeEvent &d1,
	const DssdMergeEvent &d2,
	const DssdMergeEvent &d3,
	bool *hole,
	T0Event &t0,
	TH1F *offset_window
) {
	// initialize t0 event
	t0.num = 0;
	// fill all d1 events
	for (int i = 0; i < d1.hit; ++i) {
		// if (d1.time_flag[i] != 0) continue;
		// if (d1.case_tag != 300) continue;
		t0.layer[t0.num] = 1;
		t0.flag[t0.num] = 0x1;
		t0.energy[t0.num][0] = d1.energy[i];
		t0.time[t0.num][0] = d1.time[i];
		t0.x[t0.num][0] = d1.x[i];
		t0.y[t0.num][0] = d1.y[i];
		t0.z[t0.num][0] = d1.z[i];
		t0.dssd_flag[t0.num][0] = d1.flag[i];
		++t0.num;
	}
	// fill d2 event in two cases:
	// 1. fill with d1 event if under the offset tolerance
	// 2. fill to empty slot if none match d1 event found
	for (int i = 0; i < d2.hit; ++i) {
		// if (d2.time_flag[i] != 0) continue;
		// if (d2.case_tag != 300) continue;
		bool fill = false;
		for (int j = 0; j < t0.num; ++j) {
			// jump if a d2 event has filled
			if ((t0.flag[j] & 0x2) != 0) continue;
			if ((t0.flag[j] & 0x1) == 0) continue;
			// position offset
			double xoffset = d2.x[i] - t0.x[j][0];
			double yoffset = d2.y[i] - t0.y[j][0];
			// fill to d1d2 total window
			offset_window[0].Fill(xoffset);
			offset_window[1].Fill(yoffset);
			if (fabs(xoffset) < 3.0 && fabs(yoffset) < 3.0) {
				// the angle is acceptable, fill into t0 event
				t0.layer[j]++;
				t0.flag[j] |= 0x2;
				t0.energy[j][1] = d2.energy[i];
				t0.time[j][1] = d2.time[i];
				t0.x[j][1] = d2.x[i];
				t0.y[j][1] = d2.y[i];
				t0.z[j][1] = d2.z[i];
				t0.dssd_flag[j][1] = d2.flag[i];
				fill = true;
				t0.hole[j] = hole[i];
				// this d2 event has matched, goto next one
				break;
			}
		}
		if (!fill) {
			// none of the d1 events match this d2 event,
			// put the d2 event to empty slot
			t0.layer[t0.num] = 1;
			t0.flag[t0.num] = 0x2;
			t0.energy[t0.num][0] = 0.0;
			t0.energy[t0.num][1] = d2.energy[i];
			t0.time[t0.num][0] = -1e5;
			t0.time[t0.num][1] = d2.time[i];
			t0.x[t0.num][1] = d2.x[i];
			t0.y[t0.num][1] = d2.y[i];
			t0.z[t0.num][1] = d2.z[i];
			t0.dssd_flag[t0.num][1] = d2.flag[i];
			t0.hole[t0.num] = hole[i];
			++t0.num;
		}
	}
	// fill d3 event in three cases:
	// 1. fill with d1 and d2 events under offset tolerance
	// 2. fill with only d1 event under offset tolerance
	// 3. fill with only d2 event under offset tolerance
	for (int i = 0; i < d3.hit; ++i) {
		// if (d3.time_flag[i] != 0) continue;
		// if (d3.case_tag != 300) continue;
		for (int j = 0; j < t0.num; ++j) {
			double d1d3_xoffset = d3.x[i] - t0.x[j][0];
			double d1d3_yoffset = d3.y[i] - t0.y[j][0];
			double d2d3_xoffset = d3.x[i] - t0.x[j][1];
			double d2d3_yoffset = d3.y[i] - t0.y[j][1];
			double d1d2d3_xoffset = t0.x[j][1] * 2.0 - t0.x[j][0] - d3.x[i];
			double d1d2d3_yoffset = t0.y[j][1] * 2.0 - t0.y[j][0] - d3.y[i];
			// fill d1d3 offset to total window
			offset_window[2].Fill(d1d3_xoffset);
			offset_window[3].Fill(d1d3_yoffset);
			// fill d2d3 offset to total window
			offset_window[4].Fill(d2d3_xoffset);
			offset_window[5].Fill(d2d3_yoffset);
			// fill d1d2d3 offset to total window
			offset_window[6].Fill(d1d2d3_xoffset);
			offset_window[7].Fill(d1d2d3_yoffset);
			if (t0.flag[j] == 0x3) {
				// fill d1d3 offset to flag7 window
				offset_window[10].Fill(d1d3_xoffset);
				offset_window[11].Fill(d1d3_yoffset);
				// fill d2d3 offset to flag 7 window
				offset_window[14].Fill(d2d3_xoffset);
				offset_window[15].Fill(d2d3_yoffset);
				// this slot has d1 and d2 events
				if (
					(
						fabs(d2d3_xoffset) < 4.2
						|| fabs(d1d2d3_xoffset) < 4.2
					)
					&&
					(
						fabs(d2d3_yoffset) < 4.2
						|| fabs(d1d2d3_yoffset) < 4.2
					)
				) {
					// angle check pass, fill d3 event
					t0.layer[j]++;
					t0.flag[j] |= 0x4;
					t0.energy[j][2] = d3.energy[i];
					t0.time[j][2] = d3.time[i];
					t0.x[j][2] = d3.x[i];
					t0.y[j][2] = d3.y[i];
					t0.z[j][2] = d3.z[i];
					t0.dssd_flag[j][2] = d3.flag[i];
					break;
				}
			} else if (t0.flag[j] == 0x1) {
				// fill d1d3 offset to flag5 window
				offset_window[8].Fill(d1d3_xoffset);
				offset_window[9].Fill(d1d3_yoffset);
				// this slot only contains d1 event
				if (fabs(d1d3_xoffset) < 7.0 && fabs(d1d3_yoffset) < 7.0) {
					// angle check pass, fill d3 event
					t0.layer[j]++;
					t0.flag[j] |= 0x4;
					t0.energy[j][2] = d3.energy[i];
					t0.time[j][2] = d3.time[i];
					t0.x[j][2] = d3.x[i];
					t0.y[j][2] = d3.y[i];
					t0.z[j][2] = d3.z[i];
					t0.dssd_flag[j][2] = d3.flag[i];
					break;
				}
			} else if (t0.flag[j] == 0x2) {
				// fill d2d3 offset to flag6 window
				offset_window[12].Fill(d2d3_xoffset);
				offset_window[13].Fill(d2d3_yoffset);
				// this slot only has d2 event
				if (fabs(d2d3_xoffset) < 3.0 && fabs(d2d3_yoffset) < 3.0) {
					// angle check pass, fill d3 event
					t0.layer[j]++;
					t0.flag[j] |= 0x4;
					t0.energy[j][2] = d3.energy[i];
					t0.time[j][2] = d3.time[i];
					t0.x[j][2] = d3.x[i];
					t0.y[j][2] = d3.y[i];
					t0.z[j][2] = d3.z[i];
					t0.dssd_flag[j][2] = d3.flag[i];
					break;
				}
			}
		}
	}
	// discard particles only contains the second layer
	for (int i = 0; i < t0.num; ++i) {
		while (t0.flag[i] == 0x2 && i < t0.num) {
			// this particle only contains d2 event,
			// discard it and move the following layer ahead
			for (int j = i; j < t0.num-1; ++j) {
				t0.layer[j] = t0.layer[j+1];
				t0.flag[j] = t0.flag[j+1];
				for (int k = 0; k < 3; ++k) {
					t0.energy[j][k] = t0.energy[j+1][k];
					t0.time[j][k] = t0.time[j+1][k];
					t0.x[j][k] = t0.x[j+1][k];
					t0.y[j][k] = t0.y[j+1][k];
					t0.z[j][k] = t0.z[j+1][k];
					t0.dssd_flag[j][k] = t0.dssd_flag[j+1][k];
					t0.hole[j] = t0.hole[j+1];
				}
			}
			t0.num--;
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


int T0::Track() {
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
	bool hole[8];
	d2_event.SetupInput(ipt, "d2.");
	ipt->SetBranchAddress("d2.hole", hole);
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

	// histogram of offset window between DSSD
	TH1F hist_offset[]{
		TH1F("hx12", "total x offset window of d1d2", 100, -5, 5),
		TH1F("hy12", "total y offset window of d1d2", 100, -5, 5),
		TH1F("hx13", "total x offset window of d1d3", 100, -10, 10),
		TH1F("hy13", "total y offset window of d1d3", 100, -10, 10),
		TH1F("hx23", "total x offset window of d2d3", 100, -10, 10),
		TH1F("hy23", "total y offset window of d2d3", 100, -10, 10),
		TH1F("hx123", "total x offset window of d1d2d3", 100, -10, 10),
		TH1F("hy123", "total y offset window of d1d2d3", 100, -10, 10),
		TH1F("hx13a", "total x offset window of d1d3 flag 5", 100, -10, 10),
		TH1F("hy13a", "total y offset window of d1d3 flag 5", 100, -10, 10),
		TH1F("hx13b", "total x offset window of d1d3 flag 7", 100, -10, 10),
		TH1F("hy13b", "total y offset window of d1d3 flag 7", 100, -10, 10),
		TH1F("hx23a", "total x offset window of d2d3 flag 6", 100, -5, 5),
		TH1F("hy23a", "total y offset window of d2d3 flag 6", 100, -5, 5),
		TH1F("hx23b", "total x offset window of d2d3 flag 7", 100, -5, 5),
		TH1F("hy23b", "total y offset window of d2d3 flag 7", 100, -5, 5)
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
			d1_event, d2_event, d3_event, hole,
			t0_event,
			hist_offset
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
			for (int i = 0; i < t0_event.num; ++i) {
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
	for (auto &hist : hist_offset) {
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


class Slice {
public:
	unsigned short charge;
	unsigned short mass;
	bool tail;
	int index[2];

	Slice(
		unsigned short c,
		unsigned short m,
		bool t,
		int i1,
		int i2
	)
	: charge(c)
	, mass(m)
	, tail(t)
	{
		index[0] = i1;
		index[1] = i2;
	}

	bool Connectable(const Slice &next) {
		return (
			charge == next.charge
			&& (mass == next.mass || mass == 0)
			&& index[1] == next.index[0]
			&& tail
		);
// std::cout << "Connect check:\n"
// 	<< "This layer: Z = " << charge << ", A = " << mass << ", index = " << index[1] << ", tail = " << tail << "\n"
// 	<< "Next layer: Z = " << next.charge << ", A = " << next.mass << ", index = " << next.index[0] << "\n"
// 	<< "Result: " << result << "\n";

// 		return result;
	}
};


/// @brief build DSSD-DSSD slice, e.g. T0D1-D2, T0D2-D3
/// @param[in] layer1 first layer events
/// @param[in] layer2 second layer events
/// @param[in] cuts PID cut for particle stopped in second layer
/// @param[in] tail_cuts PID cut for particle pass the second layer
/// @param[in] offset durable position offset for two layers
/// @param[out] slices built slices
///
void BuildDssdSlice(
	const DssdMergeEvent &layer1,
	const DssdMergeEvent &layer2,
	const std::vector<ParticleCut> &cuts,
	const std::vector<ParticleCut> &tail_cuts,
	double offset,
	std::vector<Slice> &slices
) {
	for (int i = 0; i < layer1.hit; ++i) {
		for (int j = 0; j < layer2.hit; ++j) {
			double xoffset = layer1.x[i] - layer2.x[j];
			double yoffset = layer1.y[i] - layer2.y[j];
			// check x, y offset
			if (fabs(xoffset) < offset && fabs(yoffset) < offset) {
				// check penatrate PID curv
				for (const auto &cut : cuts) {
					if (cut.cut->IsInside(layer2.energy[j], layer1.energy[i])) {
						// build slice
						slices.emplace_back(
							cut.charge,
							cut.mass,
							false,
							i, j
						);
					}
				}
				// check pass PID curv
				for (const auto &cut : tail_cuts) {
					if (cut.cut->IsInside(layer2.energy[j], layer1.energy[i])) {
						// build slice
						slices.emplace_back(
							cut.charge,
							cut.mass,
							true,
							i, j
						);
					}
				}
			}
		}
	}
}


/// @brief build DSSD-DSSD slice with hole, e.g. T0D1-D2, T0D2-D3
/// @param[in] layer1 first layer events
/// @param[in] layer2 second layer events
/// @param[in] cuts PID cut for particle stopped in second layer
/// @param[in] tail_cuts PID cut for particle pass the second layer
/// @param[in] hole_cuts PID cut for hole area and particle stopped
/// @param[in] hole_tail_cuts PID cut for hole area and particle passed
/// @param[in] hole1 true the first is hole
/// @param[in] hole2 true the second is hole
/// @param[in] offset durable position offset for two layers
/// @param[out] slices built slices
///
void BuildDssdSliceHole(
	const DssdMergeEvent &layer1,
	const DssdMergeEvent &layer2,
	const std::vector<ParticleCut> &cuts,
	const std::vector<ParticleCut> &tail_cuts,
	const std::vector<ParticleCut> &hole_cuts,
	const std::vector<ParticleCut> &hole_tail_cuts,
	bool *hole1,
	bool *hole2,
	double offset,
	std::vector<Slice> &slices
) {
	for (int i = 0; i < layer1.hit; ++i) {
		for (int j = 0; j < layer2.hit; ++j) {
			double xoffset = layer1.x[i] - layer2.x[j];
			double yoffset = layer1.y[i] - layer2.y[j];
			// check x, y offset
			if (fabs(xoffset) < offset && fabs(yoffset) < offset) {
				const std::vector<ParticleCut> *used_cuts = nullptr;
				const std::vector<ParticleCut> *used_tail_cuts = nullptr;
				if ((hole1 && hole1[i]) || (hole2 && hole2[j])) {
					used_cuts = &hole_cuts;
					used_tail_cuts = &hole_tail_cuts;
				} else if ((hole1 && !hole1[i]) || (hole2 && !hole2[j])) {
					used_cuts = &cuts;
					used_tail_cuts = &tail_cuts;
				} else {
					std::cerr << "Error: hole1 and hole2 invalid.\n";
					return;
				}
				// check penatrate PID curv
				for (const auto &cut : *(used_cuts)) {
					if (cut.cut->IsInside(
						layer2.energy[j], layer1.energy[i]
					)) {
						// build slice
						slices.emplace_back(
							cut.charge,
							cut.mass,
							false,
							i, j
						);
					}
				}
				// check pass PID curv
				for (const auto &cut : *(used_tail_cuts)) {
					if (cut.cut->IsInside(
						layer2.energy[j], layer1.energy[i]
					)) {
						// build slice
						slices.emplace_back(
							cut.charge,
							cut.mass,
							true,
							i, j
						);
					}
				}
			}
		}
	}
}


/// @brief build DSSD-SSD slice, e.g. T0D3-S1
/// @param[in] layer1 first layer events
/// @param[in] layer2 second layer events
/// @param[in] cuts PID cut for particle stopped in second layer
/// @param[in] tail_cuts PID cut for particle pass the second layer
/// @param[out] slices built slices
///
void BuildDssdSsdSlice(
	const DssdMergeEvent &layer1,
	const SsdEvent &layer2,
	const std::vector<ParticleCut> &cuts,
	const std::vector<ParticleCut> &tail_cuts,
	std::vector<Slice> &slices
) {
	// check valid, although energy is 0 if is invalid
	if (layer2.time < -9e4) return;
	for (int i = 0; i < layer1.hit; ++i) {
		for (const auto &cut : cuts) {
			if (cut.cut->IsInside(layer2.energy, layer1.energy[i])) {
				// build slice
				slices.emplace_back(
					cut.charge,
					cut.mass,
					false,
					i, 0
				);
			}
		}
		for (const auto &cut : tail_cuts) {
			if (cut.cut->IsInside(layer2.energy, layer1.energy[i])) {
				// build slice
				slices.emplace_back(
					cut.charge,
					cut.mass,
					true,
					i, 0
				);
			}
		}
	}
	return;
}


/// @brief build SSD-SSD slice, e.g. T0S1-S2, T0S2-S3
/// @param[in] layer1 first layer events
/// @param[in] layer2 second layer events
/// @param[in] cuts PID cut for particle stopped in second layer
/// @param[in] tail_cuts PID cut for particle pass the second layer
/// @param[out] slices built slices
///
void BuildSsdSlice(
	const SsdEvent &layer1,
	const SsdEvent &layer2,
	const std::vector<ParticleCut> &cuts,
	const std::vector<ParticleCut> &tail_cuts,
	std::vector<Slice> &slices
) {
	// check valid, although energy is 0 if is invalid
	if (layer1.time < -9e4 || layer2.time < -9e4) return;
	for (const auto &cut : cuts) {
		if (cut.cut->IsInside(layer2.energy, layer1.energy)) {
			slices.emplace_back(
				cut.charge, cut.mass,
				false,
				0, 0
			);
		}
	}
	for (const auto &cut : tail_cuts) {
		if (cut.cut->IsInside(layer2.energy, layer1.energy)) {
			slices.emplace_back(
				cut.charge, cut.mass,
				true,
				0, 0
			);
		}
	}
}


/// @brief build slices for all layers
/// @param[in] d1 T0D1 merged event
/// @param[in] d2 T0D2 merged event
/// @param[in] d3 T0D3 merged event
/// @param[in] s1 T0S1 event
/// @param[in] s2 T0S2 event
/// @param[in] s3 T0S3 event
/// @param[in] d1d2_cuts PID cut of T0D1-D2, stopped in D2
/// @param[in] d1d2_tail_cuts PID cut of T0D1-D2, pass D2
/// @param[in] d2d3_cuts PID cut of T0D2-D3, stopped in D3
/// @param[in] d2d3_tail_cuts PID cut of T0D2-D3, pass D3
/// @param[in] d3s1_cuts PID cut of T0D3-S1, stopped in S1
/// @param[in] d3s1_tail_cuts PID cut of T0D3-S1, pass S1
/// @param[in] s1s2_cuts PID cut of T0S1-S2, stopped in S2
/// @param[in] s1s2_tail_cuts PID cut of T0S1-S2, pass S2
/// @param[in] s2s3_cuts PID cut of T0S2-S3, stopped in S3
/// @param[in] s2s3_tail_cuts PID cut of T0S2-S3, pass S3
/// @param[out] slices built slices
///
void BuildSlice(
	const DssdMergeEvent &d1,
	const DssdMergeEvent &d2,
	const DssdMergeEvent &d3,
	const SsdEvent &s1,
	const SsdEvent &s2,
	const SsdEvent &s3,
	const std::vector<ParticleCut> &d1d2_cuts,
	const std::vector<ParticleCut> &d1d2_tail_cuts,
	const std::vector<ParticleCut> &d2d3_cuts,
	const std::vector<ParticleCut> &d2d3_tail_cuts,
	const std::vector<ParticleCut> &d3s1_cuts,
	const std::vector<ParticleCut> &d3s1_tail_cuts,
	const std::vector<ParticleCut> &s1s2_cuts,
	const std::vector<ParticleCut> &s1s2_tail_cuts,
	const std::vector<ParticleCut> &s2s3_cuts,
	const std::vector<ParticleCut> &s2s3_tail_cuts,
	const std::vector<ParticleCut> &d1d2_hole_cuts,
	const std::vector<ParticleCut> &d1d2_hole_tail_cuts,
	const std::vector<ParticleCut> &d2d3_hole_cuts,
	const std::vector<ParticleCut> &d2d3_hole_tail_cuts,
	bool *hole,
	std::vector<Slice> *slices
) {
	// BuildDssdSlice(d1, d2, d1d2_cuts, d1d2_tail_cuts, 4.0, slices[0]);
	// BuildDssdSlice(d2, d3, d2d3_cuts, d2d3_tail_cuts, 4.2, slices[1]);
	BuildDssdSliceHole(
		d1, d2,
		d1d2_cuts, d1d2_tail_cuts, d1d2_hole_cuts, d1d2_hole_tail_cuts,
		nullptr, hole,
		4.0,
		slices[0]
	);
	BuildDssdSliceHole(
		d2, d3,
		d2d3_cuts, d2d3_tail_cuts, d2d3_hole_cuts, d2d3_hole_tail_cuts,
		hole, nullptr,
		4.2,
		slices[1]
	);
	BuildDssdSsdSlice(d3, s1, d3s1_cuts, d3s1_tail_cuts, slices[2]);
	BuildSsdSlice(s1, s2, s1s2_cuts, s1s2_tail_cuts, slices[3]);
	BuildSsdSlice(s2, s3, s2s3_cuts, s2s3_tail_cuts, slices[4]);
}


struct SliceNode {
	// layer 0-T0D2, 1-T0D3, 2-T0S1, 3-T0S2, 4-T0S3, 5-CsI
	int layer;
	// index in slices
	int index;
	// parent node, 0 for root
	int parent;
	// flag of merge event, includes all layers' information
	// prevent to use events repeatedly
	unsigned int flag;
	// particle type get from cut, mass*10+charge, e.g. 4He is 42, 10Be is 104
	int type;
	// true if it's last layer is not in tail
	bool leaf;
};


/// @brief convert slices into tree structue
/// @param[in] slices slice nodes in two dimensional array
/// @param[in] nodes slice nodes store in vector, but connected by pointer like tree
///
void BuildSliceTree(
	std::vector<Slice> *slices,
	std::vector<SliceNode> &nodes
) {
	nodes.clear();
	// fill root node
	nodes.push_back(SliceNode{-1, -1, 0, 0, 0, false});
	// start index of each layer in nodes
	int start[6] = {
		1, 1, 1, 1, 1, 1
	};
	// offset of each layer in flag
	unsigned int offsets[6] = {
		0, 8, 16, 24, 25, 26
	};
	// fill first layer
	for (int i = 0; i < int(slices[0].size()); ++i) {
		// T0D1 flag, with offset 0
		unsigned int flag = 0x1 << slices[0][i].index[0];
		// T0D2 flag, with offset 8
		flag |= 0x100 << slices[0][i].index[1];
		// particle type
		int type = slices[0][i].mass * 10 + slices[0][i].charge;
		// insert new node
		nodes.push_back(SliceNode{
			0, i, 0, flag, type, !(slices[0][i].tail)
		});
		// increase start index
		++start[1];
	}


	for (int layer = 1; layer < 5; ++layer) {
		// intialize the next layer start value
		start[layer+1] = start[layer];
		// loop slice in this layer
		for (int i = 0; i < int(slices[layer].size()); ++i) {
			Slice &current_layer_slice = slices[layer][i];
			// search for connectable slice in nodes
			for (int j = start[layer-1]; j < start[layer]; ++j) {
// std::cout << "layer " << layer << ", slice index " << i << ", node index " << j << "\n";
				const SliceNode &node = nodes[j];
				Slice &last_layer_slice = slices[layer-1][node.index];
				if (last_layer_slice.Connectable(current_layer_slice)) {
					// flag of current layer
					unsigned int flag = 0x1 << (
						current_layer_slice.index[1] + offsets[layer+1]
					);
					// or with the previous layers
					flag |= node.flag;
					// particle type
					int type = current_layer_slice.mass*10
						+ current_layer_slice.charge;
					if (layer == 4 && current_layer_slice.tail) {
						// insert node that stop in CsI
						nodes.push_back(SliceNode{
							5, i, j, flag, type, true
						});
					} else {
						// insert new node
						nodes.push_back(SliceNode{
							layer, i, j, flag, type, !(current_layer_slice.tail)
						});
					}
					// increase start index
					++start[layer+1];
				}
			}
		}
	}

	// std::cout << "Node:\n";
	// for (size_t i = 0; i < nodes.size(); ++i) {
	// 	std::cout << i << " " << nodes[i].layer << ", "
	// 		<< nodes[i].index << ", " << nodes[i].parent
	// 		<< ", " << std::hex << nodes[i].flag << ", "
	// 		<< std::dec << nodes[i].leaf << "\n";
	// }

	return;
}


struct LeaveGroup {
	int prev;
	int index;
	int count;
};


/// @brief check whether flag is contraditory to the existing group
/// @param[in] leaves all leave nodes
/// @param[in] groups all current groups
/// @param[in] group_index index of group to check
/// @param[in] leave_flag flag of leaf to checks
/// @returns true if flag is valid, false if it's contraditory
///
bool CheckSliceFlag(
	const std::vector<SliceNode> &leaves,
	const std::vector<LeaveGroup> &groups,
	int group_index,
	unsigned int leaf_flag
) {
	for (int i = group_index; i != 0; i = groups[i].prev) {
		if ((leaves[groups[i].index].flag & leaf_flag) != 0) {
			return false;
		}
	}
	return true;
}


/// @brief Choose a group of leaves(particles)
/// @param[in] nodes all slices
/// @param[out] group selected group
/// @param[out] behe_found found 10Be and 4He in this group
/// @param[out] used_slices number of slices used in this group
///
void PickSliceGroup(
	const std::vector<SliceNode> &nodes,
	std::vector<SliceNode> &group,
	bool &behe_found,
	int &used_slices
) {
	// leaves
	std::vector<SliceNode> leaves;
	// get leaves
	for (const auto &node : nodes) {
		if (node.leaf) leaves.push_back(node);
	}

	// groups
	std::vector<LeaveGroup> groups;
	// input the root slot
	groups.push_back(LeaveGroup{-1, -1, 0});

	// insert the first leave in groups
	for (int i = 0; i < int(leaves.size()); ++i) {
		groups.push_back(LeaveGroup{0, i, 1});
	}

	// start index to add more leaves
	int start = 1;
	// tail index to add more leaves, exclusive
	int tail = groups.size();
	// combine the leaves into groups
	// if start is larger than or equal to tail, no new node has been add
	// to the groups, so, terminate the loop
	while (start < tail) {
		for (int i = start; i < tail; ++i) {
			for (int j = groups[i].index+1; j < int(leaves.size()); ++j) {
				if (CheckSliceFlag(leaves, groups, i, leaves[j].flag)) {
					// add new node to group
					groups.push_back(LeaveGroup{i, j, groups[i].count+1});
				}
			}
		}
		start = tail;
		tail = groups.size();
	}

	// // search for the group with maximum count
	// // max count
	// int max_count = 0;
	// // index with maximum count
	// int found_index = 0;
	// // loop to search
	// for (int i = 0; i < int(groups.size()); ++i) {
	// 	if (groups[i].count > max_count) {
	// 		max_count = groups[i].count;
	// 		found_index = i;
	// 	}
	// }

	behe_found = false;
	used_slices = 0;
	int found_index = 0;
	for (int i = 1; i < int(groups.size()); ++i) {
		int behe_flag = 0;
		int slice_num = 0;
		for (auto g = groups[i]; g.index != -1; g = groups[g.prev]) {
			const auto &leave = leaves[g.index];
			if (leave.type == 104) behe_flag |= 1;
			else if (leave.type == 42) behe_flag |= 2;
			slice_num += leave.layer+1;
// std::cout << "-----------------------------\n"
// 	<< "group " << i << ": " << g.index << ", prev " << g.prev << "\n"
// 	<< "type " << leave.type << ", flag " << behe_flag << "\n"
// 	<< "slice number " << slice_num << "\n";
		}
		if (behe_flag == 3) {
			if (!behe_found) {
				behe_found = true;
				used_slices = slice_num;
				found_index = i;
			} else if (slice_num > used_slices) {
				used_slices = slice_num;
				found_index = i;
			}
		} else {
			if (!behe_found && slice_num > used_slices) {
				used_slices = slice_num;
				found_index = i;
			}
		}
// std::cout << "found index " << found_index << ", behe " << (behe_found ? "Y" : "N")
// 	<< ", used slice " << used_slices << "\n";
	}

	group.clear();
	for (int i = found_index; i != 0; i = groups[i].prev) {
		group.push_back(leaves[groups[i].index]);
	}
}


int T0::SliceTrack(int supplementary) {
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
	bool hole[8];
	d2_event.SetupInput(ipt, "d2.");
	ipt->SetBranchAddress("d2.hole", hole);
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

	// supplementary T0D1 event
	DssdMergeEvent d1_suppl_event;
	// supplementary T0D2 event
	DssdMergeEvent d2_suppl_event;
	if (supplementary == 1) {
		// supplementary tag
		std::string suppl_tag = tag_.empty() ? "s1-" : tag_+"-s1-";

		// supplementary T0D1 file name
		TString d1_suppl_file_name = TString::Format(
			"%s%st0d1-merge-%s%04u.root",
			kGenerateDataPath,
			kMergeDir,
			suppl_tag.c_str(),
			run_
		);
		// add friend
		ipt->AddFriend("d1s1=tree", d1_suppl_file_name);
		// setup input branches
		d1_suppl_event.SetupInput(ipt, "d1s1.");

		// supplementary T0D2 file name
		TString d2_suppl_file_name = TString::Format(
			"%s%st0d2-merge-%s%04u.root",
			kGenerateDataPath,
			kMergeDir,
			suppl_tag.c_str(),
			run_
		);
		// add friend
		ipt->AddFriend("d2s1=tree", d2_suppl_file_name);
		// setup input branches
		d2_suppl_event.SetupInput(ipt, "d2s1.");
	}

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

	// T0D1-D2 cuts
	std::vector<ParticleCut> d1d2_cuts;
	// d1d2_cuts.push_back({1, 1, ReadCut("t0-d1d2-p1i1-1H")});
	// d1d2_cuts.push_back({1, 2, ReadCut("t0-d1d2-p1i1-2H")});
	// d1d2_cuts.push_back({1, 3, ReadCut("t0-d1d2-p1i1-3H")});
	d1d2_cuts.push_back({2, 3, ReadCut("t0-d1d2-p1i1-3He")});
	d1d2_cuts.push_back({2, 4, ReadCut("t0-d1d2-p1i1-4He")});
	d1d2_cuts.push_back({2, 6, ReadCut("t0-d1d2-p1i1-6He")});
	d1d2_cuts.push_back({3, 6, ReadCut("t0-d1d2-p1i1-6Li")});
	d1d2_cuts.push_back({3, 7, ReadCut("t0-d1d2-p1i1-7Li")});
	d1d2_cuts.push_back({3, 8, ReadCut("t0-d1d2-p1i1-8Li")});
	d1d2_cuts.push_back({3, 9, ReadCut("t0-d1d2-p1i1-9Li")});
	d1d2_cuts.push_back({4, 7, ReadCut("t0-d1d2-p1i1-7Be")});
	d1d2_cuts.push_back({4, 9, ReadCut("t0-d1d2-p1i1-9Be")});
	d1d2_cuts.push_back({4, 10, ReadCut("t0-d1d2-p1i1-10Be")});
	d1d2_cuts.push_back({5, 10, ReadCut("t0-d1d2-p1i1-10B")});
	d1d2_cuts.push_back({5, 11, ReadCut("t0-d1d2-p1i1-11B")});
	d1d2_cuts.push_back({5, 12, ReadCut("t0-d1d2-p1i1-12B")});
	d1d2_cuts.push_back({5, 13, ReadCut("t0-d1d2-p1i1-13B")});
	d1d2_cuts.push_back({6, 12, ReadCut("t0-d1d2-p1i1-12C")});
	d1d2_cuts.push_back({6, 13, ReadCut("t0-d1d2-p1i1-13C")});
	d1d2_cuts.push_back({6, 14, ReadCut("t0-d1d2-p1i1-14C")});
	// d1d2_cuts.push_back({6, 15, ReadCut("t0-d1d2-p1i1-15C")});
	// T0D1-D2 tail cuts
	std::vector<ParticleCut> d1d2_tails;
	d1d2_tails.push_back({2, 0, ReadCut("t0-d1d2-tail-p1i1-He")});
	d1d2_tails.push_back({3, 0, ReadCut("t0-d1d2-tail-p1i1-Li")});
	d1d2_tails.push_back({4, 0, ReadCut("t0-d1d2-tail-p1i1-Be")});

	// T0D2-D3 cuts
	std::vector<ParticleCut> d2d3_cuts;
	d2d3_cuts.push_back({2, 3, ReadCut("t0-d2d3-p1i1-3He")});
	d2d3_cuts.push_back({2, 4, ReadCut("t0-d2d3-p1i1-4He")});
	d2d3_cuts.push_back({2, 6, ReadCut("t0-d2d3-p1i1-6He")});
	d2d3_cuts.push_back({3, 6, ReadCut("t0-d2d3-p1i1-6Li")});
	d2d3_cuts.push_back({3, 7, ReadCut("t0-d2d3-p1i1-7Li")});
	d2d3_cuts.push_back({3, 8, ReadCut("t0-d2d3-p1i1-8Li")});
	d2d3_cuts.push_back({3, 9, ReadCut("t0-d2d3-p1i1-9Li")});
	d2d3_cuts.push_back({4, 9, ReadCut("t0-d2d3-p1i1-9Be")});
	d2d3_cuts.push_back({4, 10, ReadCut("t0-d2d3-p1i1-10Be")});
	// T0D2-D3 tail cuts
	std::vector<ParticleCut> d2d3_tails;
	d2d3_tails.push_back({2, 0, ReadCut("t0-d2d3-tail-p1i1-He")});
	d2d3_tails.push_back({3, 0, ReadCut("t0-d2d3-tail-p1i1-Li")});

	// T0D3-S1 cuts
	std::vector<ParticleCut> d3s1_cuts;
	d3s1_cuts.push_back({2, 3, ReadCut("t0-d3s1-p1i1-3He")});
	d3s1_cuts.push_back({2, 4, ReadCut("t0-d3s1-p1i1-4He")});
	d3s1_cuts.push_back({2, 6, ReadCut("t0-d3s1-p1i1-6He")});
	d3s1_cuts.push_back({3, 6, ReadCut("t0-d3s1-p1i1-6Li")});
	d3s1_cuts.push_back({3, 7, ReadCut("t0-d3s1-p1i1-7Li")});
	d3s1_cuts.push_back({3, 8, ReadCut("t0-d3s1-p1i1-8Li")});
	d3s1_cuts.push_back({3, 9, ReadCut("t0-d3s1-p1i1-9Li")});
	// T0D3-S1 tail cuts
	std::vector<ParticleCut> d3s1_tails;
	d3s1_tails.push_back({2, 0, ReadCut("t0-d3s1-tail-p1i1-He")});

	// T0S1-S2 cuts
	std::vector<ParticleCut> s1s2_cuts;
	s1s2_cuts.push_back({2, 3, ReadCut("t0-s1s2-p1i1-3He")});
	s1s2_cuts.push_back({2, 4, ReadCut("t0-s1s2-p1i1-4He")});
	s1s2_cuts.push_back({2, 6, ReadCut("t0-s1s2-p1i1-6He")});
	// T0S1-S2 cuts
	std::vector<ParticleCut> s1s2_tails;
	s1s2_tails.push_back({2, 0, ReadCut("t0-s1s2-tail-p1i1-He")});
	s1s2_tails.push_back({2, 6, ReadCut("t0-s1s2-tail-p1i1-6He")});

	// T0S2-S3 cuts
	std::vector<ParticleCut> s2s3_cuts;
	s2s3_cuts.push_back({2, 3, ReadCut("t0-s2s3-p1i1-3He")});
	s2s3_cuts.push_back({2, 4, ReadCut("t0-s2s3-p1i1-4He")});
	s2s3_cuts.push_back({2, 6, ReadCut("t0-s2s3-p1i1-6He")});
	// T0S2-S3 pass cuts
	std::vector<ParticleCut> s2s3_tails;
	s2s3_tails.push_back({2, 4, ReadCut("t0-s2s3-tail-p1i1-4He")});
	s2s3_tails.push_back({2, 6, ReadCut("t0-s2s3-tail-p1i1-6He")});

	// T0D1-D2 hole cuts
	std::vector<ParticleCut> d1d2_hole_cuts;
	d1d2_hole_cuts.push_back({2, 4, ReadCut("t0-d1d2-hole-p1i1-4He")});
	d1d2_hole_cuts.push_back({4, 10, ReadCut("t0-d1d2-hole-p1i1-10Be")});
	// T0D1-D2 hole tail cuts
	std::vector<ParticleCut> d1d2_hole_tail_cuts;
	d1d2_hole_tail_cuts.push_back({2, 0, ReadCut("t0-d1d2-hole-tail-p1i1-He")});
	d1d2_hole_tail_cuts.push_back({4, 0, ReadCut("t0-d1d2-hole-tail-p1i1-Be")});

	// T0D2-D3 hole cuts
	std::vector<ParticleCut> d2d3_hole_cuts;
	d2d3_hole_cuts.push_back({2, 4, ReadCut("t0-d2d3-hole-p1i1-4He")});
	d2d3_hole_cuts.push_back({4, 10, ReadCut("t0-d2d3-hole-p1i1-10Be")});
	// T0D2-D3 hole tail cuts
	std::vector<ParticleCut> d2d3_hole_tail_cuts;
	d2d3_hole_tail_cuts.push_back({2, 0, ReadCut("t0-d2d3-hole-tail-p1i1-He")});

	// for statistics
	long long normal_count = 0;
	long long d1_suppl_count = 0;
	long long d2_suppl_count = 0;
	long long d1d2_suppl_count = 0;


	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Tracking t0   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
	// for (long long entry : test_entries) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);
// std::cout << "============= entry " << entry << " ==============\n";


		// picked slices
		std::vector<Slice> slices[5];
		// picked nodes
		std::vector<SliceNode> nodes;
		// picked nodess group
		std::vector<SliceNode> group;
		// used events
		DssdMergeEvent *used_d1_event = &d1_event;
		DssdMergeEvent *used_d2_event = &d2_event;

		// slices without supplementary
		std::vector<Slice> normal_slices[5];
		// nodes without supplementary
		std::vector<SliceNode> normal_nodes;
		// temporary picked nodes group
		std::vector<SliceNode> normal_group;
		// found 10Be and 4He in normal group ?
		bool normal_found_behe;
		// number of used slices in normal group
		int normal_used_slices;
		// build slice from DSSD and SSD events
		BuildSlice(
			d1_event, d2_event, d3_event,
			s1_event, s2_event, s3_event,
			d1d2_cuts, d1d2_tails,
			d2d3_cuts, d2d3_tails,
			d3s1_cuts, d3s1_tails,
			s1s2_cuts, s1s2_tails,
			s2s3_cuts, s2s3_tails,
			d1d2_hole_cuts, d1d2_hole_tail_cuts,
			d2d3_hole_cuts, d2d3_hole_tail_cuts,
			hole,
			normal_slices
		);
		// convert to tree structure
		BuildSliceTree(normal_slices, normal_nodes);
		// pick one group
		PickSliceGroup(
			normal_nodes, normal_group,
			normal_found_behe, normal_used_slices
		);

		if (supplementary == 0) {
			for (int i = 0; i < 5; ++i) slices[i] = normal_slices[i];
			nodes = normal_nodes;
			group = normal_group;
			used_d1_event = &d1_event;
			used_d2_event = &d2_event;
		} else if (supplementary == 1) {
			// slices with d1 and d2 supplementary events
			std::vector<Slice> d1d2_suppl_slices[5];
			// nodes with d1 and d2 supplementary events
			std::vector<SliceNode> d1d2_suppl_nodes;
			// picked nodes group with d1 and d2 supplementary events
			std::vector<SliceNode> d1d2_suppl_group;
			// found 10Be and 4He in d1 and d2 supplementary group ?
			bool d1d2_suppl_found_behe;
			// number of used slices in and d2 supplementary group
			int d1d2_suppl_used_slices;
			// build slices with d1 and d2 supplementary events
			BuildSlice(
				d1_suppl_event, d2_suppl_event, d3_event,
				s1_event, s2_event, s3_event,
				d1d2_cuts, d1d2_tails,
				d2d3_cuts, d2d3_tails,
				d3s1_cuts, d3s1_tails,
				s1s2_cuts, s1s2_tails,
				s2s3_cuts, s2s3_tails,
				d1d2_hole_cuts, d1d2_hole_tail_cuts,
				d2d3_hole_cuts, d2d3_hole_tail_cuts,
				hole,
				d1d2_suppl_slices
			);
			// convert to tree structure
			BuildSliceTree(d1d2_suppl_slices, d1d2_suppl_nodes);
			// pick one group
			PickSliceGroup(
				d1d2_suppl_nodes, d1d2_suppl_group,
				d1d2_suppl_found_behe, d1d2_suppl_used_slices
			);

			// evaluate different group and pick the bset one
			// selected flag, 0-normal, 1-d1_suppl, 2-d2_suppl, 3-d1d2_suppl
			int selected = 0;

			if (normal_found_behe) {
				selected = 0;
				++normal_count;
			// } else if (d2_suppl_found_behe) {
			// 	selected = 2;
			// 	++d2_suppl_count;
			} else if (d1d2_suppl_found_behe) {
				selected = 3;
				++d1d2_suppl_count;
			// } else if (d1_suppl_found_behe) {
			// 	selected = 1;
			// 	++d1_suppl_count;
			} else {
				size_t max_count = normal_group.size();
				int max_slices = normal_used_slices;

				if (
					d1d2_suppl_group.size() > max_count
					|| (
						d1d2_suppl_group.size() == max_count
						&& d1d2_suppl_used_slices > max_slices
					)
				) {
					selected = 3;
					max_count = d1d2_suppl_group.size();
					max_slices = d1d2_suppl_used_slices;
				}
			}
			if (selected == 0) {
				for (int i = 0; i < 5; ++i) slices[i] = normal_slices[i];
				nodes = normal_nodes;
				group = normal_group;
				used_d1_event = &d1_event;
				used_d2_event = &d2_event;
				++normal_count;
			} else {
				for (int i = 0; i < 5; ++i) slices[i] = d1d2_suppl_slices[i];
				nodes = d1d2_suppl_nodes;
				group = d1d2_suppl_group;
				used_d1_event = &d1_suppl_event;
				used_d2_event = &d2_suppl_event;
			}
		}


// std::cout << "-------------------------\nEntry: " << entry << "\n";
// std::cout << "Slices:\n";
// for (int i = 0; i < 5; ++i) {
// 	for (const auto &slice : slices[i]) {
// 		std::cout << "Slice layer " << i
// 			<< ": Z = " << slice.charge << ", A = " << slice.mass
// 			<< ", tail = " << slice.tail << ", index = " << slice.index[0]
// 			<< ", " << slice.index[1] << "\n";
// 	}
// }
// std::cout << "Slice tree:\n";
// for (const auto &node : nodes) {
// 	std::cout << "layer " << node.layer << ", index = " << node.index
// 		<< ", parent " << node.parent << ", flag = " << std::hex << node.flag
// 		<< std::dec << ", type = " << node.type << ", leaf = "
// 		<< node.leaf << "\n";
// }
// std::cout << "Slice group:\n";
// for (const auto &node : group) {
// 	std::cout << "layer " << node.layer << ", index = " << node.index
// 		<< ", parent " << node.parent << ", flag = " << std::hex << node.flag
// 		<< std::dec << ", type = " << node.type << ", leaf = "
// 		<< node.leaf << "\n";
// }

		// fill output event
		t0_event.num = int(group.size());
		// if (t0_event.num > 8) t0_event.num = 8;
		for (int i = 0; i < t0_event.num; ++i) {
			int detect_layer = group[i].layer == 5 ? 4 : group[i].layer;
			if (group[i].index >= int(slices[detect_layer].size())) {
				std::cerr
					<< "Error: group[i] index over range in run " << run_
					<< ", entry " << entry
					<< ", layer " << detect_layer
					<< ", index is " << group[i].index
					<< ", range " << slices[detect_layer].size()
					<< " .\n";
				return -1;
			}
			Slice *slice = &(slices[detect_layer][group[i].index]);

			t0_event.layer[i] = group[i].layer + 1;
			t0_event.flag[i] = group[i].layer == 0 ? 0x3 : 0x7;
			t0_event.ssd_flag = group[i].flag >> 24;
			t0_event.charge[i] = slice->charge;
			t0_event.mass[i] = slice->mass;

			SliceNode *node = &(group[i]);
			for (int layer = detect_layer; layer >= 0; --layer) {
				if (node->index >= int(slices[layer].size())) {
					std::cerr
						<< "Error: node index over range in run " << run_
						<< ", entry " << entry
						<< ", layer " << layer
						<< ", index is " << node->index
						<< ", range " << slices[layer].size()
						<< " .\n";
					return -1;
				}
				slice = &(slices[layer][node->index]);
				if (layer == 0) {
					// fill T0D1 information
					int d1_index = slice->index[0];
					if (d1_index >= used_d1_event->hit) {
						std::cerr << "Error: D1 index over hit in run "
							<< run_ << ", entry "
							<< entry << ", index " << d1_index
							<< ", hit " << used_d1_event->hit << "\n";
						return -1;
					}
					t0_event.energy[i][0] = used_d1_event->energy[d1_index];
					t0_event.time[i][0] = used_d1_event->time[d1_index];
					t0_event.x[i][0] = used_d1_event->x[d1_index];
					t0_event.y[i][0] = used_d1_event->y[d1_index];
					t0_event.z[i][0] = used_d1_event->z[d1_index];
					t0_event.dssd_flag[i][0] = used_d1_event->flag[d1_index];
					// fill T0D2 information
					int d2_index = slice->index[1];
					if (d2_index >= used_d2_event->hit) {
						std::cerr << "Error: D2 index over hit in run "
							<< run_ << ", entry "
							<< entry << ", index " << d2_index
							<< ", hit " << used_d2_event->hit << "\n";
						return -1;
					}
					t0_event.energy[i][1] = used_d2_event->energy[d2_index];
					t0_event.time[i][1] = used_d2_event->time[d2_index];
					t0_event.x[i][1] = used_d2_event->x[d2_index];
					t0_event.y[i][1] = used_d2_event->y[d2_index];
					t0_event.z[i][1] = used_d2_event->z[d2_index];
					t0_event.dssd_flag[i][1] = used_d2_event->flag[d2_index];
					t0_event.hole[i] = hole[d2_index];
				} else if (layer == 1) {
					// fill T0D3 information
					int d3_index = slice->index[1];
					if (d3_index >= d3_event.hit) {
						std::cerr << "Error: D1 index over hit in run "
							<< run_ << ", entry "
							<< entry << ", index " << d3_index
							<< ", hit " << d3_event.hit << " .\n";
						return -1;
					}
					t0_event.energy[i][2] = d3_event.energy[d3_index];
					t0_event.time[i][2] = d3_event.time[d3_index];
					t0_event.x[i][2] = d3_event.x[d3_index];
					t0_event.y[i][2] = d3_event.y[d3_index];
					t0_event.z[i][2] = d3_event.z[d3_index];
					t0_event.dssd_flag[i][2] = d3_event.flag[d3_index];
				} else if (layer == 2) {
					// fill T0S1 information
					t0_event.ssd_energy[0] = s1_event.energy;
					t0_event.ssd_time[0] = s1_event.time;
				} else if (layer == 3) {
					t0_event.ssd_energy[1] = s2_event.energy;
					t0_event.ssd_time[1] = s2_event.time;
				} else if (layer == 4) {
					t0_event.ssd_energy[2] = s3_event.energy;
					t0_event.ssd_time[2] = s3_event.time;
				}
				node = &(nodes[node->parent]);
			}

			// Hole correct: Some events still in normal cut even though
			// they're in hole pixel, maybe I can remark them as not hole
			// if (t0_event.hole[i]) {
			// 	if (t0_event.charge[i] == 4 && t0_event.mass[i] == 10) {
			// 		// 10Be hole event
			// 		if (t0_event.layer[i] == 1) {
			// 			// stop in T0D2
			// 			if (d1d2_cuts[8].cut->IsInside(
			// 				t0_event.energy[i][1], t0_event.energy[i][0]
			// 			)) {
			// 				t0_event.mass[i] = 9;
			// 				t0_event.hole[i] = false;
			// 			} else if (d1d2_cuts[9].cut->IsInside(
			// 				t0_event.energy[i][1], t0_event.energy[i][0]
			// 			)) {
			// 				t0_event.mass[i] = 10;
			// 				t0_event.hole[i] = false;
			// 			}
			// 		} else if (t0_event.layer[i] == 2) {
			// 			// stop in T0D3
			// 			if (
			// 				d1d2_tails[2].cut->IsInside(
			// 					t0_event.energy[i][1], t0_event.energy[i][0]
			// 				)
			// 				&&
			// 				d2d3_cuts[7].cut->IsInside(
			// 					t0_event.energy[i][2], t0_event.energy[i][1]
			// 				)
			// 			) {
			// 				t0_event.mass[i] = 9;
			// 				t0_event.hole[i] = false;
			// 			} else if (
			// 				d1d2_tails[2].cut->IsInside(
			// 					t0_event.energy[i][1], t0_event.energy[i][0]
			// 				)
			// 				&&
			// 				d2d3_cuts[8].cut->IsInside(
			// 					t0_event.energy[i][2], t0_event.energy[i][1]
			// 				)
			// 			) {
			// 				t0_event.mass[i] = 10;
			// 				t0_event.hole[i] = false;
			// 			}
			// 		}

			// 	} else if (t0_event.charge[i] == 2 && t0_event.mass[i] == 4) {
			// 		// 4He hole event
			// 		if (t0_event.layer[i] == 1) {
			// 			// stop in T0D2
			// 			if (d1d2_cuts[1].cut->IsInside(
			// 				t0_event.energy[i][1], t0_event.energy[i][0]
			// 			)) {
			// 				t0_event.hole[i] = false;
			// 			}
			// 		} else if (t0_event.layer[i] == 2) {
			// 			// stop in T0D3
			// 			if (
			// 				d1d2_tails[0].cut->IsInside(
			// 					t0_event.energy[i][1], t0_event.energy[i][0]
			// 				)
			// 				&&
			// 				d2d3_cuts[1].cut->IsInside(
			// 					t0_event.energy[i][2], t0_event.energy[i][1]
			// 				)
			// 			) {
			// 				t0_event.hole[i] = false;
			// 			}
			// 		} else if (t0_event.layer[i] > 2) {
			// 			// stop in SSD or CsI
			// 			if (
			// 				d1d2_tails[0].cut->IsInside(
			// 					t0_event.energy[i][1], t0_event.energy[i][0]
			// 				)
			// 				&&
			// 				d2d3_tails[0].cut->IsInside(
			// 					t0_event.energy[i][2], t0_event.energy[i][1]
			// 				)
			// 			) {
			// 				t0_event.hole[i] = false;
			// 			}
			// 		}
			// 	}
			// }
		}

		opt.Fill();
	}

	// show finish
	printf("\b\b\b\b100%%\n");

	std::vector<bool> identify;
	identify.push_back(true);

	// save tree
	opt.Write();
	t0_file.WriteObject(&identify, "identify");
	// close files
	t0_file.Close();
	d1_file.Close();

	if (supplementary == 1) {
		std::cout << "Normal " << normal_count << "\n"
			<< "D1 supplementary " << d1_suppl_count << "\n"
			<< "D2 supplementary " << d2_suppl_count << "\n"
			<< "D1D2 supplementary " << d1d2_suppl_count << "\n";
	}

	return 0;
}


int DssdParticleIdentify(
	const T0Event &t0,
	size_t index,
	const std::vector<ParticleCut> &d1d2_cuts,
	const std::vector<ParticleCut> &d1d2_tails,
	const std::vector<ParticleCut> &d1d2_center_cuts,
	const std::vector<ParticleCut> &d1d2_center_tails,
	const std::vector<ParticleCut> &d2d3_cuts,
	const std::vector<ParticleCut> &d2d3_center_cuts,
	ParticleTypeEvent &type
) {
	bool center = false;
	// if (
	// 	t0.x[index][1] > -9 && t0.x[index][1] < -3
	// 	&& t0.y[index][1] > -2 && t0.y[index][1] < 7
	// ) {
	// 	center = true;
	// }
	if (t0.flag[index] == 0x3) {
		const std::vector<ParticleCut> &cuts =
			center ? d1d2_center_cuts : d1d2_cuts;
		// particle pass 2 DSSD
		for (const auto &cut : cuts) {
			if (cut.cut->IsInside(t0.energy[index][1], t0.energy[index][0])) {
				// particle is in cuts, recored the charge and mass number
				type.charge[index] = cut.charge;
				type.mass[index] = cut.mass;
				type.layer[index] = 1;
				return 1;
			}
		}
	} else if (t0.flag[index] == 0x7) {
		const std::vector<ParticleCut> &cuts =
			center ? d2d3_center_cuts : d2d3_cuts;
		const std::vector<ParticleCut> &tails =
			center ? d1d2_center_tails : d1d2_tails;
		// particle pass all 3 DSSD
		for (size_t i = 0; i < cuts.size(); ++i) {
			if (
				cuts[i].cut->IsInside(
					t0.energy[index][2], t0.energy[index][1]
				)
				&& tails[i].cut->IsInside(
					t0.energy[index][1], t0.energy[index][0]
				)
			) {
				// particle is in cuts, recored the charge and mass number
				type.charge[index] = cuts[i].charge;
				type.mass[index] = cuts[i].mass;
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
	const std::vector<ParticleCut> &d2d3_tails,
	const std::vector<ParticleCut> &d2d3_center_tails,
	const std::vector<ParticleCut> &d3s1_cuts,
	const std::vector<ParticleCut> &d3s1_tails,
	const std::vector<ParticleCut> &s1s2_cuts,
	const std::vector<ParticleCut> &s1s2_tails,
	const std::vector<ParticleCut> &s2s3_cuts,
	const std::vector<ParticleCut> &s2s3_tails,
	const std::unique_ptr<TCutG> &s1s2_cross_cut,
	ParticleTypeEvent &type
) {
	bool center = false;
	// if (
	// 	t0.x[index][1] > -9 && t0.x[index][1] < -3
	// 	&& t0.y[index][1] > -2 && t0.y[index][1] < 7
	// ) {
	// 	center = true;
	// }

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
		const std::vector<ParticleCut> &tails =
			center ? d2d3_center_tails : d2d3_tails;
		// particle stop in the first SSD
		for (size_t i = 0; i < d3s1_cuts.size(); ++i) {
			if (
				d3s1_cuts[i].cut->IsInside(
					t0.ssd_energy[0], t0.energy[index][2]
				)
				&& tails[i].cut->IsInside(
					t0.energy[index][2], t0.energy[index][1]
				)
			) {
				// particle is in cuts, record the charge and mass
				type.charge[index] = d3s1_cuts[i].charge;
				type.mass[index] = d3s1_cuts[i].mass;
				type.layer[index] = 3;
				return 3;
			}
		}
	} else if (ssd_flag == 0x3) {
		// particle stop in the second SSD
		for (size_t i = 0; i < s1s2_cuts.size(); ++i) {
			if (
				s1s2_cuts[i].cut->IsInside(
					t0.ssd_energy[1], t0.ssd_energy[0]
				)
				&& d3s1_tails[i].cut->IsInside(
					t0.ssd_energy[0], t0.energy[index][2]
				)
			) {
				// particle is in cuts, record the charge and mass
				type.charge[index] = s1s2_cuts[i].charge;
				type.mass[index] = s1s2_cuts[i].mass;
				type.layer[index] = 4;
				return 4;
			}
		}
	} else if (ssd_flag == 0x7) {
		// particle stop in the third SSD or CsI(Tl)
		for (size_t i = 0; i < s2s3_cuts.size(); ++i) {
			if (
				s2s3_cuts[i].cut->IsInside(
					t0.ssd_energy[2], t0.ssd_energy[1]
				)
				&& s1s2_tails[i].cut->IsInside(
					t0.ssd_energy[1], t0.ssd_energy[0]
				)
			) {
				type.charge[index] = s2s3_cuts[i].charge;
				type.mass[index] = s2s3_cuts[i].mass;
				type.layer[index] = 5;
				return 5;
			}
		}
		for (const auto &cut : s2s3_tails) {
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
	// d1d2_cuts.push_back({1, 1, ReadCut("t0-d1d2-p1i1-1H")});
	// d1d2_cuts.push_back({1, 2, ReadCut("t0-d1d2-p1i1-2H")});
	// d1d2_cuts.push_back({1, 3, ReadCut("t0-d1d2-p1i1-3H")});
	d1d2_cuts.push_back({2, 3, ReadCut("t0-d1d2-p1i1-3He")});
	d1d2_cuts.push_back({2, 4, ReadCut("t0-d1d2-p1i1-4He")});
	d1d2_cuts.push_back({2, 6, ReadCut("t0-d1d2-p1i1-6He")});
	d1d2_cuts.push_back({3, 6, ReadCut("t0-d1d2-p1i1-6Li")});
	d1d2_cuts.push_back({3, 7, ReadCut("t0-d1d2-p1i1-7Li")});
	d1d2_cuts.push_back({4, 7, ReadCut("t0-d1d2-p1i1-7Be")});
	d1d2_cuts.push_back({4, 9, ReadCut("t0-d1d2-p1i1-9Be")});
	d1d2_cuts.push_back({4, 10, ReadCut("t0-d1d2-p1i1-10Be")});
	d1d2_cuts.push_back({5, 10, ReadCut("t0-d1d2-p1i1-10B")});
	d1d2_cuts.push_back({5, 11, ReadCut("t0-d1d2-p1i1-11B")});
	d1d2_cuts.push_back({5, 12, ReadCut("t0-d1d2-p1i1-12B")});
	d1d2_cuts.push_back({5, 13, ReadCut("t0-d1d2-p1i1-13B")});
	d1d2_cuts.push_back({6, 12, ReadCut("t0-d1d2-p1i1-12C")});
	d1d2_cuts.push_back({6, 13, ReadCut("t0-d1d2-p1i1-13C")});
	d1d2_cuts.push_back({6, 14, ReadCut("t0-d1d2-p1i1-14C")});
	// d1d2_cuts.push_back({6, 15, ReadCut("t0-d1d2-p1i1-15C")});
	// T0D1-D2 tail cuts
	std::vector<ParticleCut> d1d2_tails;
	d1d2_tails.push_back({2, 0, ReadCut("t0-d1d2-tail-p1i1-He")});
	d1d2_tails.push_back({3, 0, ReadCut("t0-d1d2-tail-p1i1-Li")});
	d1d2_tails.push_back({4, 0, ReadCut("t0-d1d2-tail-p1i1-Be")});
	// T0D1-D2 center cuts
	std::vector<ParticleCut> d1d2_center_cuts;
	// d1d2_center_cuts.push_back({2, 4, ReadCut("t0-d1d2-center-He")});
	// d1d2_center_cuts.push_back({4, 10, ReadCut("t0-d1d2-center-Be")});
	// T0D1-D2 center tail cuts
	std::vector<ParticleCut> d1d2_center_tails;
	// d1d2_center_tails.push_back({2, 0, ReadCut("t0-d1d2-tail-center-He")});
	// d1d2_center_tails.push_back({4, 0, ReadCut("t0-d1d2-tail-center-Be")});

	// T0D2-D3 cuts
	std::vector<ParticleCut> d2d3_cuts;
	d2d3_cuts.push_back({2, 4, ReadCut("t0-d2d3-p1i1-4He")});
	d2d3_cuts.push_back({3, 7, ReadCut("t0-d2d3-p1i1-7Li")});
	d2d3_cuts.push_back({4, 10, ReadCut("t0-d2d3-p1i1-10Be")});
	// T0D2-D3 tail cuts
	std::vector<ParticleCut> d2d3_tails;
	d2d3_tails.push_back({2, 0, ReadCut("t0-d2d3-tail-p1i1-He")});
	// T0D2-D3 center cuts
	std::vector<ParticleCut> d2d3_center_cuts;
	// d2d3_center_cuts.push_back({2, 4, ReadCut("t0-d2d3-center-He")});
	// d2d3_center_cuts.push_back({4, 10, ReadCut("t0-d2d3-center-Be")});
	// T0D2-D3 center tail cuts
	std::vector<ParticleCut> d2d3_center_tails;
	// d2d3_center_tails.push_back({2, 0, ReadCut("t0-d2d3-tail-center-He")});

	// T0D3-S1 cuts
	std::vector<ParticleCut> d3s1_cuts;
	d3s1_cuts.push_back({2, 4, ReadCut("t0-d3s1-p1i1-4He")});
	// T0D3-S1 tail cuts
	std::vector<ParticleCut> d3s1_tails;
	d3s1_tails.push_back({2, 0, ReadCut("t0-d3s1-tail-p1i1-He")});
	// T0S1-S2 cuts
	std::vector<ParticleCut> s1s2_cuts;
	s1s2_cuts.push_back({2, 4, ReadCut("t0-s1s2-p1i1-4He")});
	// T0S1-S2 cuts
	std::vector<ParticleCut> s1s2_tails;
	s1s2_tails.push_back({2, 0, ReadCut("t0-s1s2-tail-p1i1-He")});
	// T0S2-S3 cuts
	std::vector<ParticleCut> s2s3_cuts;
	s2s3_cuts.push_back({2, 4, ReadCut("t0-s2s3-p1i1-4He")});
	// T0S2-S3 pass cuts
	std::vector<ParticleCut> s2s3_tails;
	s2s3_tails.push_back({2, 4, ReadCut("t0-s2s3-tail-p1i1-He")});
	// special cut of S1-S2 interaction
	std::unique_ptr<TCutG> s1s2_cross_cut{ReadCut("t0-s1s2-cross")};

	//statistics
	long long total = 0;
	long long id11 = 0;
	long long id21 = 0;
	long long id22 = 0;
	long long id31 = 0;
	long long id32 = 0;
	long long id33 = 0;

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
		for (int i = 0; i < type_event.num; ++i) {
			type_event.charge[i] = 0;
			type_event.mass[i] = 0;
			type_event.layer[i] = -1;
		}
		// identify particles
		if (t0_event.num == 1) {
			++total;
			int identify = DssdParticleIdentify(
				t0_event, 0,
				d1d2_cuts, d1d2_tails, d1d2_center_cuts, d1d2_center_tails,
				d2d3_cuts, d2d3_center_cuts,
				type_event
			);
			if (
				identify <= 0
				&& t0_event.flag[0] == 0x7
				&& t0_event.ssd_flag != 0
			) {
				identify = SsdParticleIdentify(
					t0_event, 0,
					d2d3_tails, d2d3_center_tails, d3s1_cuts,
					d3s1_tails, s1s2_cuts,
					s1s2_tails, s2s3_cuts,
					s2s3_tails, s1s2_cross_cut,
					type_event
				);
			}
			if (identify > 0) ++id11;
		} else if (t0_event.num == 2) {
			++total;
			int identify0 = DssdParticleIdentify(
				t0_event, 0,
				d1d2_cuts, d1d2_tails, d1d2_center_cuts, d1d2_center_tails,
				d2d3_cuts, d2d3_center_cuts,
				type_event
			);
			int identify1 = DssdParticleIdentify(
				t0_event, 1,
				d1d2_cuts, d1d2_tails, d1d2_center_cuts, d1d2_center_tails,
				d2d3_cuts, d2d3_center_cuts,
				type_event
			);
			if (identify0 <= 0 && identify1 > 0 && t0_event.flag[0] == 0x7) {
				int identify = SsdParticleIdentify(
					t0_event, 0,
					d2d3_tails, d2d3_center_tails, d3s1_cuts,
					d3s1_tails, s1s2_cuts,
					s1s2_tails, s2s3_cuts,
					s2s3_tails, s1s2_cross_cut,
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
					d2d3_tails, d2d3_center_tails, d3s1_cuts,
					d3s1_tails, s1s2_cuts,
					s1s2_tails, s2s3_cuts,
					s2s3_tails, s1s2_cross_cut,
					type_event
				);
				if (identify > 0) ++id22;
				else ++id21;
			} else if (identify0 > 0 && identify1 > 0) {
				++id22;
			} else if (identify0 <= 0 && identify1 > 0) {
				++id21;
			} else if (identify0 > 0 && identify1 <= 0) {
				++id21;
			}
		} else if (t0_event.num == 3) {
			++total;
			int identify0 = DssdParticleIdentify(
				t0_event, 0,
				d1d2_cuts, d1d2_tails, d1d2_center_cuts, d1d2_center_tails,
				d2d3_cuts, d2d3_center_cuts,
				type_event
			);
			int identify1 = DssdParticleIdentify(
				t0_event, 1,
				d1d2_cuts, d1d2_tails, d1d2_center_cuts, d1d2_center_tails,
				d2d3_cuts, d2d3_center_cuts,
				type_event
			);
			int identify2 = DssdParticleIdentify(
				t0_event, 2,
				d1d2_cuts, d1d2_tails, d1d2_center_cuts, d1d2_center_tails,
				d2d3_cuts, d2d3_center_cuts,
				type_event
			);
			if (identify0 > 0 && identify1 > 0 && identify2 > 0) {
				++id33;
			} else if (
				identify0 <= 0 && identify1 > 0 && identify2 > 0
				&& t0_event.flag[0] == 0x7
			) {
				int identify = SsdParticleIdentify(
					t0_event, 0,
					d2d3_tails, d2d3_center_tails, d3s1_cuts,
					d3s1_tails, s1s2_cuts,
					s1s2_tails, s2s3_cuts,
					s2s3_tails, s1s2_cross_cut,
					type_event
				);
				if (identify > 0) ++id33;
				else ++id32;
			} else if (
				identify0 > 0 && identify1 <= 0 && identify2 > 0
				&& t0_event.flag[1] == 0x7
			) {
				int identify = SsdParticleIdentify(
					t0_event, 1,
					d2d3_tails, d2d3_center_tails, d3s1_cuts,
					d3s1_tails, s1s2_cuts,
					s1s2_tails, s2s3_cuts,
					s2s3_tails, s1s2_cross_cut,
					type_event
				);
				if (identify > 0) ++id33;
				else ++id32;
			} else if (
				identify0 > 0 && identify1 > 0 && identify2 <= 0
				&& t0_event.flag[2] == 0x7
			) {
				int identify = SsdParticleIdentify(
					t0_event, 2,
					d2d3_tails, d2d3_center_tails, d3s1_cuts,
					d3s1_tails, s1s2_cuts,
					s1s2_tails, s2s3_cuts,
					s2s3_tails, s1s2_cross_cut,
					type_event
				);
				if (identify > 0) ++id33;
				else ++id32;
			} else if (identify0 > 0 || identify1 > 0 || identify2 > 0) {
				++id31;
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
		<< "  " << double(id22) / double(total) << "\n"
		<< "Identify 1/3 " << id31 << " / " << total
		<< "  " << double(id31) / double(total) << "\n"
		<< "Identify 2/3 " << id32 << " / " << total
		<< "  " << double(id32) / double(total) << "\n"
		<< "Idnetify 3/3 " << id33 << " / " << total
		<< "  " << double(id33) / double(total) << "\n";
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
	// {0, 5, 6, 14, 36'000.0, 52'000.0, 0.0},
	// d2d3
	{1, 2, 2, 4, 5000.0, 9000.0, 55'000.0},
	{1, 3, 3, 7, 10'000.0, 16'000.0, 55'000.0},
	{1, 4, 4, 10, 20'000.0, 26'000.0, 55'000.0},
	// d3s1
	// {2, 1, 1, 2, 1900.0, 2400.0, 85'000.0},
	{2, 2, 2, 4, 5200.0, 9500.0, 85'000.0},
	// s1s2
	// {3, 0, 1, 1, 3500.0, 5000.0, 95'000.0},
	{3, 2, 2, 4, 13'500.0, 25'000.0, 95'000.0},
	// s2s3
	// {4, 0, 1, 1, 3500.0, 5500.0, 120'000.0},
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


class CurveFit {
public:
	CurveFit(): calculator_("t0", "4He") {}
	double operator()(double *x, double *par) const {
		double de = par[0] + par[1] * x[0];
		double e = calculator_.Energy(0, de);
		return (e - par[2]) / par[3];
	}
private:
	elc::DeltaEnergyCalculator calculator_;
};


int T0::Calibrate(unsigned int end_run) {
	// input T0 telescope chain
	TChain t0_chain("t0", "chain of T0 events");
	for (unsigned int i = run_; i <= end_run; ++i) {
		if (i == 628) continue;
		t0_chain.AddFile(TString::Format(
			"%s%st0-telescope-%s%04u.root/tree",
			kGenerateDataPath,
			kTelescopeDir,
			tag_.empty() ? "" : (tag_+"-").c_str(),
			i
		));
	}
	// input T0 particle type chain
	TChain type_chain("type", "chain of particle type events");
	for (unsigned int i = run_; i <= end_run; ++i) {
		if (i == 628) continue;
		type_chain.AddFile(TString::Format(
			"%s%st0-particle-type-%s%04u.root/tree",
			kGenerateDataPath,
			kParticleIdentifyDir,
			tag_.empty() ? "" : (tag_+"-").c_str(),
			i
		));
	}
	// add friend
	t0_chain.AddFriend(&type_chain, "type");
	// input T0 telescope event
	T0Event t0_event;
	// input T0 particle type event
	ParticleTypeEvent type_event;
	// setup input branches
	t0_event.SetupInput(&t0_chain);
	type_event.SetupInput(&t0_chain, "type.");

	// output calibration root file name
	TString calibration_file_name;
	calibration_file_name.Form(
		"%s%st0-calibration-%s%04u-%04u.root",
		// "%s%st0-calibration-%s%04u.root",
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


	// total number of entries
	long long entries = t0_chain.GetEntries();
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
		t0_chain.GetEntry(entry);
		for (int i = 0; i < t0_event.num; ++i) {
			if (type_event.layer[i] == 1) {
				double de = t0_event.energy[i][0];
				double e = t0_event.energy[i][1];
				if (
					type_event.charge[i] == 2
					&& type_event.mass[i] == 4
				) {
					if (de > 4000.0 && de < 8000.0) {
						e_vs_de_offset.AddPoint(de, e);
					}
				} else if (
					type_event.charge[i] == 3
					&& type_event.mass[i] == 7
				) {
					if (de > 8000.0 && de < 13'000.0) {
						e_vs_de_offset.AddPoint(de, e);
					}
				} else if (
					type_event.charge[i] == 4
					&& type_event.mass[i] == 10
				) {
					if (de > 13'000.0 && de < 24'000.0) {
						e_vs_de_offset.AddPoint(de, e);
					}
				}
			} else if (type_event.layer[i] == 2) {
				double de = t0_event.energy[i][1];
				double e = t0_event.energy[i][2];
				if (
					type_event.charge[i] == 2
					&& type_event.mass[i] == 4
				) {
					if (de > 5000.0 && de < 9000.0) {
						e_vs_de_offset.AddPoint(de+55'000.0, e);
					}
				} else if (
					type_event.charge[i] == 3
					&& type_event.mass[i] == 7
				) {
					if (de > 10'000.0 && de < 16'000.0) {
						e_vs_de_offset.AddPoint(de+55'000.0, e);
					}
				} else if (
					type_event.charge[i] == 4
					&& type_event.mass[i] == 10
				) {
					if (de > 20'000.0 && de < 26'000.0) {
						e_vs_de_offset.AddPoint(de+55'000.0, e);
					}
				}
			} else if (type_event.layer[i] == 3) {
				double de = t0_event.energy[i][2];
				double e = t0_event.ssd_energy[0];
				if (
					type_event.charge[i] == 2
					&& type_event.mass[i] == 4
				) {
					if (de > 5200.0 && de < 9500.0) {
						e_vs_de_offset.AddPoint(de+85'000.0, e);
					}
				}
			} else if (type_event.layer[i] == 4) {
				double de = t0_event.ssd_energy[0];
				double e = t0_event.ssd_energy[1];
				if (
					type_event.charge[i] == 2
					&& type_event.mass[i] == 4
				) {
					if (de > 13'500.0 && de < 25'000.0) {
						e_vs_de_offset.AddPoint(de+95'000.0, e);
					}
				}
			} else if (type_event.layer[i] == 5) {
				double de = t0_event.ssd_energy[1];
				double e = t0_event.ssd_energy[2];
				if (de > 13'500.0 && de < 22'500.0) {
					e_vs_de_offset.AddPoint(de+120'000.0, e);
				}
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
		return -1;
	}

	// save graph
	calibration_file.cd();
	e_vs_de_offset.Write("gcali");
	// close files
	calibration_file.Close();

	return 0;
}



int T0::ReadCalibrateParameters(unsigned int end_run) {
	// parameters file name
	TString file_name;
	if (end_run == 9999) {
		file_name.Form(
			"%s%s%s-calibration-param.txt",
			kGenerateDataPath,
			kCalibrationDir,
			name_.c_str()
		);
	} else if (end_run == run_) {
		file_name.Form(
			"%s%s%s-calibration-param-%04u.txt",
			kGenerateDataPath,
			kCalibrationDir,
			name_.c_str(),
			run_
		);
	} else {
		file_name.Form(
			"%s%s%s-calibration-param-%04u-%04u.txt",
			kGenerateDataPath,
			kCalibrationDir,
			name_.c_str(),
			run_,
			end_run
		);
	}
	// parameters file
	std::ifstream fin(file_name.Data());
	// check file
	if (!fin.good()) {
		std::cerr << "Error: Read calibrate parameters from "
			<< file_name << " failed.\n";
		return -1;
	}
	// read parameters
	for (size_t i = 0; i < Layers(); ++i) {
		fin >> cali_params_[i*2] >> cali_params_[i*2+1];
	}
	// close file
	fin.close();
	return 0;
}


int T0::WriteCalibrateParameters(unsigned int end_run) const {
	// parameters file name
	TString file_name;
	if (end_run == 9999) {
		file_name.Form(
			"%s%s%s-calibration-param.txt",
			kGenerateDataPath,
			kCalibrationDir,
			name_.c_str()
		);
	} else if (end_run == run_) {
		file_name.Form(
			"%s%s%s-calibration-param-%04u.txt",
			kGenerateDataPath,
			kCalibrationDir,
			name_.c_str(),
			run_
		);
	} else {
		file_name.Form(
			"%s%s%s-calibration-param-%04u-%04u.txt",
			kGenerateDataPath,
			kCalibrationDir,
			name_.c_str(),
			run_,
			end_run
		);
	}
	// parameters file
	std::ofstream fout(file_name.Data());
	// check file
	if (!fout.good()) {
		std::cerr << "Error: Open calibrate file "
			<< file_name << " failed.\n";
		return -1;
	}
	// write parameters
	for (size_t i = 0; i < Layers(); ++i) {
		fout << cali_params_[i*2] << " " << cali_params_[i*2+1] << "\n";
	}
	// close file
	fout.close();
	return 0;
}


void T0::CalibrateResult(T0Event &t0) {
	for (int i = 0; i < t0.num; ++i) {
		for (int j = 0; j < 3; ++j) {
			t0.energy[i][j] = CaliEnergy(j, t0.energy[i][j]);
		}
	}
	for (int i = 0; i < 3; ++i) {
		t0.ssd_energy[i] = CaliEnergy(i+3, t0.ssd_energy[i]);
	}
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
	// input telescope event
	T0Event t0_event;
	// setup input branches
	t0_event.SetupInput(ipt);

	// has identify?
	std::vector<bool> *identify =
		(std::vector<bool>*)telescope_file.Get("identify");

	// input type event
	ParticleTypeEvent type_event;
	if (!identify || !(identify->at(0))) {
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
		// setup input branches
		type_event.SetupInput(ipt, "type.");
	} else {
		ipt->SetBranchAddress("charge", type_event.charge);
		ipt->SetBranchAddress("mass", type_event.mass);
		ipt->SetBranchAddress("layer", type_event.layer);
	}

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
	bool hole[8];
	bool has_4he;
	bool has_10be;
	// setup output branches
	particle_event.SetupOutput(&opt);
	opt.Branch("hole", hole, "hole[num]/O");
	opt.Branch("has_4he", &has_4he, "he/O");
	opt.Branch("has_10be", &has_10be, "be/O");

	// read calibrate parameters
	if (ReadCalibrateParameters() && ReadCalibrateParameters(9999)) {
		std::cerr << "Errro: Read calibrate parameters failed.\n";
		return -1;
	}

	// CsI energy calculator
	elc::CsiEnergyCalculator csi_calculator("4He");
	// delta energy calculator
	elc::DeltaEnergyCalculator delta_4He_calculator("t0", "4He");
	elc::DeltaEnergyCalculator delta_10Be_calculator("t0", "10Be");

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
		has_4he = has_10be = false;
		for (int i = 0; i < t0_event.num; ++i) {
			// jump confuesd particles
			if (type_event.charge[i] <= 0 || type_event.mass[i] <= 0) continue;
			if (type_event.layer[i] < 1 || type_event.layer[i] > 6) continue;

			// fill charge number
			particle_event.charge[particle_event.num] = type_event.charge[i];
			// fill mass number
			particle_event.mass[particle_event.num] = type_event.mass[i];
			// fill layer
			particle_event.layer[particle_event.num] = type_event.layer[i];

			double energy = TotalEnergy(
				t0_event, type_event, i, csi_calculator, delta_10Be_calculator
			);
			// fill energy
			particle_event.energy[particle_event.num] = energy;
			particle_event.time[particle_event.num] = t0_event.time[i][0];

			// set particle position, T0D1 may has the best resolution
			particle_event.x[particle_event.num] = t0_event.x[i][0];
			particle_event.y[particle_event.num] = t0_event.y[i][0];
			particle_event.z[particle_event.num] = t0_event.z[i][0];

			// leave momentum empty
			particle_event.px[particle_event.num] = 0.0;
			particle_event.py[particle_event.num] = 0.0;
			particle_event.pz[particle_event.num] = 0.0;

			// other information
			particle_event.status[particle_event.num] = t0_event.status[i];
			particle_event.index[particle_event.num] = i;
			hole[particle_event.num] = t0_event.hole[i];
			if (type_event.mass[i] == 10 && type_event.charge[i] == 4) {
				has_10be = true;
			}
			if (type_event.mass[i] == 4 && type_event.charge[i] == 2) {
				has_4he = true;
			}
			particle_event.straight[particle_event.num] = t0_event.straight[i];

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
	const elc::CsiEnergyCalculator &csi_calculator,
	const elc::DeltaEnergyCalculator &/*delta_calculator*/
) const {
	// total energy
	double energy = CaliEnergy(0, t0.energy[index][0]);
	// DSSD layers particle hit
	if (type.layer[index] >= 1) {
		// if (type.charge[index] == 4 && type.mass[index] == 10) {
			// energy += delta_calculator.Energy(0, energy);
		// } else {
			energy += CaliEnergy(1, t0.energy[index][1]);
		// }
	}
	if (type.layer[index] >= 2) energy += CaliEnergy(2, t0.energy[index][2]);
	// return energy if particle stop at T0D2
	if (type.layer[index] == 1) return energy;

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


int T0::ShowCalibration() {
	// telescope file name
	TString tele_file_name = TString::Format(
		"%s%st0-telescope-%s%04u.root",
		kGenerateDataPath,
		kTelescopeDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// input file
	TFile tele_file(tele_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)tele_file.Get("tree");
	// particle type file name
	TString pid_file_name = TString::Format(
		"%s%st0-particle-type-%s%04u.root",
		kGenerateDataPath,
		kParticleIdentifyDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// add friend
	ipt->AddFriend("pid=tree", pid_file_name);
	// input event
	T0Event t0_event;
	ParticleTypeEvent type_event;
	// setup input branches
	t0_event.SetupInput(ipt);
	type_event.SetupInput(ipt, "pid.");

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%st0-cali-%s%04u.root",
		kGenerateDataPath,
		kShowDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// 1D histogram of T0D1 energy difference, all identified particles
	TH1F d1_identified_de("d1ide", "energy difference", 1000, -10, 10);
	// 1D histogram of T0D1 energy difference, only 4He
	TH1F d1_he_de("d1hede", "energy difference", 1000, -10, 10);
	// 1D histogram of T0D2 energy difference, all identified particles
	TH1F d2_identified_de("d2ide", "energy difference", 1000, -10, 10);
	// 1D histogram of T0D2 energy difference, only 4He
	TH1F d2_he_de("d2hede", "energy difference", 1000, -10, 10);

	// read calibration parameters from file
	if (ReadCalibrateParameters(run_)) {
		std::cerr << "Error: Read calibration parameters failed.\n";
		return -1;
	}

	// range-energy calculator
	std::vector<elc::RangeEnergyCalculator> calculators;
	calculators.emplace_back("4He", "Si");
	calculators.emplace_back("6Li", "Si");
	calculators.emplace_back("7Li", "Si");
	calculators.emplace_back("9Be", "Si");
	calculators.emplace_back("10Be", "Si");

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Generating calibrate result   0%%");
	fflush(stdout);
	// loop to calculate calibrated energy
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		ipt->GetEntry(entry);
		// calculate DSSD energy
		for (int i = 0; i < t0_event.num; ++i) {
			for (size_t j = 0; j < 3; ++j) {
				t0_event.energy[i][j] *= cali_params_[j*2+1];
				t0_event.energy[i][j] += cali_params_[j*2];
			}
		}
		// calculate SSD energy
		for (size_t j = 0; j < 3; ++j) {
			t0_event.ssd_energy[j] *= cali_params_[(j+3)*2+1];
			t0_event.ssd_energy[j] += cali_params_[(j+3)*2];
		}
		// check flag
		for (int i = 0; i < t0_event.num; ++i) {
			if (t0_event.flag[i] != 0x3) continue;
			if (type_event.mass[i] == 0 || type_event.charge[i] == 0) continue;
			size_t index = 0;
			// sum of energy
			double total_energy = t0_event.energy[i][0] + t0_event.energy[i][1];
			if (type_event.mass[i] == 4) {
				index = 0;
			} else if (type_event.mass[i] == 6) {
				index = 1;
			} else if (type_event.mass[i] == 7 && type_event.charge[i] == 3) {
				index = 2;
			} else if (type_event.mass[i] == 9 && type_event.charge[i] == 4) {
				index = 3;
			} else if (type_event.mass[i] == 10 && type_event.charge[i] == 4) {
				index = 4;
			}
			// incident range
			double range = calculators[index].Range(total_energy);
			// d2 range
			double d2_range = range - t0_thickness[0];
			// d2 energy
			double d2_energy = calculators[index].Energy(d2_range);
			// d1 energy
			double d1_energy = total_energy - d2_energy;
			// fill to histogram
			d1_identified_de.Fill(d1_energy - t0_event.energy[i][0]);
			d2_identified_de.Fill(d2_energy - t0_event.energy[i][1]);
			if (index == 0) {
				d1_he_de.Fill(d1_energy - t0_event.energy[i][0]);
				d2_he_de.Fill(d2_energy - t0_event.energy[i][1]);
			}

		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save histograms
	opf.cd();
	d1_he_de.Write();
	d2_he_de.Write();
	d1_identified_de.Write();
	d2_identified_de.Write();
	// close files
	opf.Close();
	tele_file.Close();
	return 0;
}


constexpr double comb_energy_threshold = 1000;
constexpr double single_comb_time_threshold = 200;
constexpr double comb_time_threshold = 400;

struct StripGroup {
	double energy;
	double time;
	double strip;
	unsigned int flag;
	int points;
	int status;
};

void SearchStripCombinations(
	DssdFundamentalEvent &event,
	TH1F &fadt,
	TH1F &badt,
	std::vector<StripGroup> &front_strips,
	std::vector<StripGroup> &back_strips
) {
	// for convenience
	int &fhit = event.front_hit;
	int &bhit = event.back_hit;
	double *fe = event.front_energy;
	double *be = event.back_energy;
	double *ft = event.front_time;
	double *bt = event.back_time;
	unsigned short *fs = event.front_strip;
	unsigned short *bs = event.back_strip;

	// search for front strips
	for (int i = 0; i < fhit; ++i) {
		front_strips.push_back({fe[i], ft[i], double(fs[i]), 1u<<i, 0, 0});
		for (int j = i+1; j < fhit; ++j) {
			if (abs(fs[i]-fs[j]) == 1) {
				fadt.Fill(ft[i]-ft[j]);
			}
			if (
				abs(fs[i]-fs[j]) == 1
				&& fabs(ft[i]-ft[j]) < single_comb_time_threshold
			) {
				front_strips.push_back({
					fe[i]+fe[j],
					ft[i],
					double(fs[i])+fe[j]/(fe[i]+fe[j])*(fs[j]-fs[i]),
					(1u << i) | (1u << j),
					-1,
					1
				});
			}
		}
	}
	// search for back strips
	for (int i = 0; i < bhit; ++i) {
		back_strips.push_back({be[i], bt[i], double(bs[i]), 0x100u<<i, 0, 0});
		for (int j = i+1; j < bhit; ++j) {
			if (abs(bs[i]-bs[j]) == 1) {
				badt.Fill(bt[i]-bt[j]);
			}
			if (
				abs(bs[i]-bs[j]) == 1
				&& fabs(bt[i]-bt[j]) < single_comb_time_threshold
			) {
				back_strips.push_back({
					be[i]+be[j],
					bt[i],
					double(bs[i])+be[j]/(be[i]+be[j])*(bs[j]-bs[i]),
					(0x100u << i) | (0x100u << j),
					-1,
					1
				});
			}
		}
	}

// std::cout << "\n-----------------strip combinations------------------------\n";
// std::cout << std::hex;
// for (const auto &comb : front_strips) {
// 	std::cout << "status " << comb.status << " flag " << comb.flag << "\n"
// 		<< "energy " << comb.energy << " time " << comb.time << "\n";
// }
// for (const auto &comb : back_strips) {
// 	std::cout << "status " << comb.status << " flag " << comb.flag << "\n"
// 		<< "energy " << comb.energy << " time " << comb.time << "\n";
// }
}



ROOT::Math::XYZVector CalculatePosition(size_t index, double fs, double bs) {
	if (index == 0) {
		// T0D1
		double x = bs + 0.5 - 32;
		double y = fs + 0.5 - 32;
		ROOT::Math::XYZVector result (x, y, 100.0);
		return result;
	} else if (index == 1) {
		// T0D2
		double x = 2.0 * (fs + 0.5) - 32.0 - 1.04;
		double y = 2.0 * (bs + 0.5) - 32.0 - 0.90;
		ROOT::Math::XYZVector result(x, y, 111.76);
		return result;
	} else if (index == 2) {
		// T0D3
		double x = 2.0 * (fs + 0.5) - 32.0 - 1.16;
		double y = 2.0 * (bs + 0.5) - 32.0 - 1.0;
		return ROOT::Math::XYZVector(x, y, 123.52);
	}
	return {0.0, 0.0, 0.0};
}


struct SideCombination {
	int num;
	double energy[2];
	double time[2];
	double x[2];
	double y[2];
	int points;
	unsigned int flag[2];
	int status[2];
};


std::vector<SideCombination> SearchSideCombinations(
	const std::vector<StripGroup> &front,
	const std::vector<StripGroup> &back,
	size_t index,
	TH1F &sde,
	TH1F &sdt
) {
	// result
	std::vector<SideCombination> result;
	// search for separating events
	for (size_t i = 0; i < front.size(); ++i) {
		for (size_t j = 0; j < back.size(); ++j) {
			if (fabs(front[i].energy - back[j].energy) > comb_energy_threshold) continue;
			if (fabs(front[i].time - back[j].time) > single_comb_time_threshold) continue;
			sde.Fill(front[i].energy - back[j].energy);
			sdt.Fill(front[i].time - back[j].time);
			int points = 0;
			if (fabs(front[i].energy - back[j].energy) < 500) {
				points += 4;
			} else if (fabs(front[i].energy - back[j].energy) < 1000) {
				points += 2;
			} else if (fabs(front[i].energy - back[j].energy) < 2000) {
				points += 1;
			}
			if (fabs(front[i].time - back[j].time) < 20) {
				points += 5;
			} else if (fabs(front[i].time - back[j].time) < 80) {
				points += 3;
			} else if (fabs(front[i].time - back[j].time) < 150) {
				points += 2;
			} else if (fabs(front[i].time - back[j].time) < 200) {
				points += 1;
			}
			auto position = CalculatePosition(
				index, front[i].strip, back[j].strip
			);
			result.push_back({
				1,
				front[i].energy, 0.0,
				front[i].time, 0.0,
				position.X(), 0.0,
				position.Y(), 0.0,
				points + front[i].points + back[j].points,
				front[i].flag | back[j].flag, 0,
				front[i].status + back[j].status*4, 0
			});
		}
	}
	size_t single_size = result.size();
	// search for f2b1 binding events
	for (size_t i = 0; i < front.size(); ++i) {
		if (front[i].status != 0) continue;
		for (size_t j = i+1; j < front.size(); ++j) {
			if (front[j].status != 0) continue;
			for (size_t k = 0; k < back.size(); ++k) {
				if (back[k].status != 0) continue;
				// check energy
				double de =
					fabs(front[i].energy + front[j].energy - back[k].energy);
				if (de > comb_energy_threshold) continue;
				// check time
				double dt1 = fabs(front[i].time - front[j].time);
				double dt2 = fabs(front[i].time - back[k].time);
				if (
					dt1 > single_comb_time_threshold
					|| dt2 > single_comb_time_threshold
				) continue;
				int points = 0;
				// evaluate energy points
				if (de < 500) points += 4;
				else if (de < 1000) points += 2;
				else if (de < 2000) points += 1;
				// evaluate time points
				if (dt1 < 20) points += 3;
				else if (dt1 < 80) points += 2;
				else if (dt1 < 150) points += 1;
				// else if (dt1 < 200) points += 1;
				if (dt2 < 20) points += 3;
				else if (dt2 < 80) points += 2;
				else if (dt2 < 150) points += 1;
				// else if (dt2 < 200) points += 1;
				auto position1 = CalculatePosition(
					index, front[i].strip, back[k].strip
				);
				auto position2 = CalculatePosition(
					index, front[j].strip, back[k].strip
				);
				result.push_back({
					2,
					front[i].energy,
					front[j].energy,
					front[i].time,
					front[j].time,
					position1.X(),
					position2.X(),
					position1.Y(),
					position2.Y(),
					points + front[i].points + back[j].points,
					front[i].flag | back[k].flag,
					front[j].flag | back[k].flag,
					2*4,
					2*4
				});
			}
		}
	}
	// search for f1b2 binding events
	for (size_t i = 0; i < front.size(); ++i) {
		if (front[i].status != 0) continue;
		for (size_t j = 0; j < back.size(); ++j) {
			if (back[j].status != 0) continue;
			for (size_t k = j+1; k < back.size(); ++k) {
				if (back[k].status != 0) continue;
				// check energy
				double de =
					fabs(front[i].energy - back[j].energy - back[k].energy);
				if (de > comb_energy_threshold) continue;
				// check time
				double dt1 = fabs(back[j].time - back[k].time);
				double dt2 = fabs(front[i].time - back[j].time);
				if (
					dt1 > single_comb_time_threshold
					|| dt2 > single_comb_time_threshold
				) continue;
				int points = 0;
				// evaluate energy points
				if (de < 500) points += 4;
				else if (de < 1000) points += 2;
				else if (de < 2000) points += 1;
				// evaluate time points
				if (dt1 < 20) points += 3;
				else if (dt1 < 80) points += 2;
				else if (dt1 < 150) points += 1;
				// else if (dt1 < 200) points += 1;
				if (dt2 < 20) points += 3;
				else if (dt2 < 80) points += 2;
				else if (dt2 < 150) points += 1;
				// else if (dt2 < 200) points += 1;
				auto position1 = CalculatePosition(
					index, front[i].strip, back[j].strip
				);
				auto position2 = CalculatePosition(
					index, front[i].strip, back[k].strip
				);
				result.push_back({
					2,
					back[j].energy,
					back[k].energy,
					back[j].time,
					back[k].time,
					position1.X(),
					position2.X(),
					position1.Y(),
					position2.Y(),
					points + front[i].points + back[j].points,
					front[i].flag | back[j].flag,
					front[i].flag | back[k].flag,
					2,
					2
				});
			}
		}
	}
	// search for f2b2 binding events
	for (size_t i = 0; i < single_size; ++i) {
		for (size_t j = i+1; j < single_size; ++j) {
			// ignore binding events
			if (result[i].num != 1 || result[j].num != 1) continue;
			// ignore adjacent strips events
			if (result[i].status[0] != 0 || result[j].status[0] != 0) continue;
			// check conflict
			if ((result[i].flag[0] & result[j].flag[0]) != 0) continue;
			// check time
			if (fabs(result[i].time[0] - result[j].time[0]) > 200) continue;
			result.push_back({
				4,
				result[i].energy[0],
				result[j].energy[0],
				result[i].time[0],
				result[j].time[0],
				result[i].x[0],
				result[j].x[0],
				result[i].y[0],
				result[j].y[0],
				(result[i].points + result[j].points) / 2,
				result[i].flag[0],
				result[j].flag[0],
				0,
				0
			});
		}
	}

// std::cout << "\nsingle size " << std::dec << single_size << "\n";
// std::cout << "---------------side combinations-----------\n" << std::hex;
// for (const auto & comb : result) {
// 	std::cout << comb.num << " " << comb.flag[0] << " " << comb.flag[1] << "\n";
// }
	return result;
}

struct Combination {
	int points;
	bool operator<(const Combination &other) const {
		return points < other.points;
	}
};

struct TrackCombination {
	double energy[2];
	double time[2];
	double x[2];
	double y[2];
	unsigned int flag[2];
	int status;
	int points;

	bool operator<(const TrackCombination& other) const {
		return points < other.points;
	}
};


struct TrackBindingCombination {
	double energy[2][2];
	double time[2][2];
	double x[2][2];
	double y[2][2];
	unsigned int flag[2][2];
	int status[2];
	int points;

	bool operator<(const TrackBindingCombination& other) const {
		return points < other.points;
	}
};


struct D3Combination {
	double energy;
	double time;
	double x;
	double y;
	unsigned int flag;
	int status;
	int points;
	int d1d2_index;

	bool operator<(const D3Combination &other) const {
		return points < other.points;
	}
};


int D1D2ParticleIdentify(
	double d1_energy,
	double d2_energy,
	const std::vector<ParticleCut> &d1d2_cuts
) {
	// particle stop in T0D2
	for (const auto &cut : d1d2_cuts) {
		if (cut.cut->IsInside(d2_energy, d1_energy)) {
			return cut.charge + cut.mass*10;
		}
	}
	return 0;
}


int SearchTrackCombinations(
	const std::vector<SideCombination>& d1_combinations,
	const std::vector<SideCombination>& d2_combinations,
	const std::vector<SideCombination>& d3_combinations,
	// const std::vector<ParticleCut> &d1d2_cuts,
	const std::vector<ParticleCut> &d1d2_tails,
	// const std::vector<ParticleCut> &d2d3_cuts,
	unsigned int d1_flag,
	unsigned int d2_flag,
	unsigned int d3_flag,
	TH1F &d1d2dt, TH1F &d1d2dx, TH1F &d1d2dy,
	TH1F &d1d3dt, TH1F &d1d3dx, TH1F &d1d3dy,
	TH1F &d2d3dt, TH1F &d2d3dx, TH1F &d2d3dy,
	TH1F &d1d2_points, TH1F &bind_points, TH1F &d3_points,
	T0Event &t0
) {
	// search for seperating events
	// queue of combinations
	std::priority_queue<TrackCombination> queue;
	for (const auto &d1 : d1_combinations) {
		if (d1.num != 1) continue;
		for (const auto &d2 : d2_combinations) {
			if (d2.num != 1) continue;
			if (fabs(d1.x[0]-d2.x[0]) > 4 || fabs(d1.y[0]-d2.y[0]) > 4) continue;
			if (fabs(d1.time[0]-d2.time[0]) > comb_time_threshold) continue;
			d1d2dt.Fill(d1.time[0] - d2.time[0]);
			d1d2dx.Fill(d1.x[0] - d2.x[0]);
			d1d2dy.Fill(d1.y[0] - d2.y[0]);

			int points = 0;

			if (fabs(d1.x[0]-d2.x[0]) < 1.0) points += 3;
			else if (fabs(d1.x[0]-d2.x[0]) < 2.0) points += 2;
			else if (fabs(d1.x[0]-d2.x[0]) < 3.0) points += 1;

			if (fabs(d1.y[0]-d2.y[0]) < 1.0) points += 3;
			else if (fabs(d1.y[0]-d2.y[0]) < 2.0) points += 2;
			else if (fabs(d1.y[0]-d2.y[0]) < 3.0) points += 1;

			if (fabs(d1.time[0] - d2.time[0]) < 20.0) points += 5;
			else if (fabs(d1.time[0] - d2.time[0]) < 80.0) points += 3;
			else if (fabs(d1.time[0] - d2.time[0]) < 150.0) points += 2;
			else if (fabs(d1.time[0] - d2.time[0]) < 200.0) points += 1;
			else points -= 1;

			TrackCombination comb{
				{d1.energy[0], d2.energy[0]},
				{d1.time[0], d2.time[0]},
				{d1.x[0], d2.x[0]},
				{d1.y[0], d2.y[0]},
				{d1.flag[0], d2.flag[0]},
				d1.status[0] + d2.status[0]*16,
				points + d1.points + d2.points
			};
			d1d2_points.Fill(comb.points);
			queue.push(comb);
		}
	}
	// initialize
	t0.num = 0;
	while (!queue.empty() && d1_flag != 0 && d2_flag != 0 && t0.num < 8) {
		TrackCombination comb = queue.top();
		queue.pop();
		if (comb.points < 20) continue;
		if ((d1_flag & comb.flag[0]) != comb.flag[0]) continue;
		if ((d2_flag & comb.flag[1]) != comb.flag[1]) continue;
// std::cout << "\n" << t0.num << " d1 flag " << std::hex << comb.flag[0] << " / " << d1_flag << "\n"
// 	<< comb.flag[1] << " / " << d2_flag << "\n";
		d1_flag &= ~comb.flag[0];
		d2_flag &= ~comb.flag[1];
		t0.layer[t0.num] = 2;
		t0.flag[t0.num] = 0x3;
		t0.energy[t0.num][0] = comb.energy[0];
		t0.energy[t0.num][1] = comb.energy[1];
		t0.time[t0.num][0] = comb.time[0];
		t0.time[t0.num][1] = comb.time[1];
		t0.x[t0.num][0] = comb.x[0];
		t0.x[t0.num][1] = comb.x[1];
		t0.y[t0.num][0] = comb.y[0];
		t0.y[t0.num][1] = comb.y[1];
		t0.z[t0.num][0] = 100.0;
		t0.z[t0.num][1] = 111.76;
		t0.status[t0.num] = comb.status;
		t0.points[t0.num] = comb.points;
		t0.dssd_flag[t0.num][0] = comb.flag[0];
		t0.dssd_flag[t0.num][1] = comb.flag[1];
		++t0.num;
	}

	// search for binding events
	std::priority_queue<TrackBindingCombination> bind_queue;
	for (const auto &d1 : d1_combinations) {
		if (d1.num != 2 && d1.num != 4) continue;
		for (const auto &d2 : d2_combinations) {
			if (d2.num != 2 && d2.num != 4) continue;
			if (d1.num == 4 && d2.num == 4) continue;
			if (fabs(d1.time[0]-d2.time[0]) > comb_time_threshold) continue;
			if (fabs(d1.time[1]-d2.time[1]) > comb_time_threshold) continue;
			if (
				!(
					fabs(d1.x[0]-d2.x[0]) < 4
					&& fabs(d1.y[0]-d2.y[0]) < 4
					&& fabs(d1.x[1]-d2.x[1]) < 4
					&& fabs(d1.y[1]-d2.y[1]) < 4
				)
				&&
				!(
					fabs(d1.x[0]-d2.x[1]) < 4
					&& fabs(d1.y[0]-d2.y[1]) < 4
					&& fabs(d1.x[1]-d2.x[0]) < 4
					&& fabs(d1.y[1]-d2.y[0]) < 4
				)
			) continue;
			d1d2dt.Fill(d1.time[0] - d2.time[0]);
			d1d2dx.Fill(d1.x[0] - d2.x[0]);
			d1d2dy.Fill(d1.y[0] - d2.y[0]);

			int points = 0;

			if (fabs(d1.time[0] - d2.time[0]) < 20.0) points += 5;
			else if (fabs(d1.time[0] - d2.time[0]) < 80.0) points += 3;
			else if (fabs(d1.time[0] - d2.time[0]) < 150.0) points += 2;
			else if (fabs(d1.time[0] - d2.time[0]) < 200.0) points += 1;
			else points -= 1;

			int position_points1 = 0;
			if (fabs(d1.x[0]-d2.x[0]) < 1.0) position_points1 += 3;
			else if (fabs(d1.x[0]-d2.x[0]) < 2.0) position_points1 += 2;
			else if (fabs(d1.x[0]-d2.x[0]) < 3.0) position_points1 += 1;

			if (fabs(d1.y[0]-d2.y[0]) < 1.0) position_points1 += 3;
			else if (fabs(d1.y[0]-d2.y[0]) < 2.0) position_points1 += 2;
			else if (fabs(d1.y[0]-d2.y[0]) < 3.0) position_points1 += 1;

			if (fabs(d1.x[1]-d2.x[1]) < 1.0) position_points1 += 3;
			else if (fabs(d1.x[1]-d2.x[1]) < 2.0) position_points1 += 2;
			else if (fabs(d1.x[1]-d2.x[1]) < 3.0) position_points1 += 1;

			if (fabs(d1.y[1]-d2.y[1]) < 1.0) position_points1 += 3;
			else if (fabs(d1.y[1]-d2.y[1]) < 2.0) position_points1 += 2;
			else if (fabs(d1.y[1]-d2.y[1]) < 3.0) position_points1 += 1;


			int position_points2 = 0;
			if (fabs(d1.x[0]-d2.x[1]) < 1.0) position_points2 += 3;
			else if (fabs(d1.x[0]-d2.x[1]) < 2.0) position_points2 += 2;
			else if (fabs(d1.x[0]-d2.x[1]) < 3.0) position_points2 += 1;

			if (fabs(d1.y[0]-d2.y[1]) < 1.0) position_points2 += 3;
			else if (fabs(d1.y[0]-d2.y[1]) < 2.0) position_points2 += 2;
			else if (fabs(d1.y[0]-d2.y[1]) < 3.0) position_points2 += 1;

			if (fabs(d1.x[1]-d2.x[0]) < 1.0) position_points2 += 3;
			else if (fabs(d1.x[1]-d2.x[0]) < 2.0) position_points2 += 2;
			else if (fabs(d1.x[1]-d2.x[0]) < 3.0) position_points2 += 1;

			if (fabs(d1.y[1]-d2.y[0]) < 1.0) position_points2 += 3;
			else if (fabs(d1.y[1]-d2.y[0]) < 2.0) position_points2 += 2;
			else if (fabs(d1.y[1]-d2.y[0]) < 3.0) position_points2 += 1;

			// combination 1
			TrackBindingCombination comb1{
				{d1.energy[0], d2.energy[0], d1.energy[1], d2.energy[1]},
				{d1.time[0], d2.time[0], d1.time[1], d2.time[1]},
				{d1.x[0], d2.x[0], d1.x[1], d2.x[1]},
				{d1.y[0], d2.y[0], d1.y[1], d2.y[1]},
				{d1.flag[0], d1.flag[1], d2.flag[0], d2.flag[1]},
				d1.status[0]+d2.status[0]*16, d1.status[1]+d2.status[1]*16,
				points + position_points1 + d1.points + d2.points
			};
			if (!(
				d1.num == 2
				&& d2.num == 2
				&& (
					d1.status[0] == d2.status[0]
					|| d1.status[1] == d2.status[1]
				)
			)) {
				bind_points.Fill(comb1.points);
				bind_queue.push(comb1);
			}
			// combination 2
			TrackBindingCombination comb2{
				{d1.energy[0], d2.energy[1], d1.energy[1], d2.energy[0]},
				{d1.time[0], d2.time[1], d1.time[1], d2.time[0]},
				{d1.x[0], d2.x[1], d1.x[1], d2.x[0]},
				{d1.y[0], d2.y[1], d1.y[1], d2.y[0]},
				{d1.flag[0], d1.flag[1], d2.flag[1], d2.flag[0]},
				d1.status[0]+d2.status[1]*16, d1.status[1]+d2.status[0]*16,
				points + position_points2 + d1.points + d2.points
			};
			if (!(
				d1.num == 2
				&& d2.num == 2
				&& (
					d1.status[0] == d2.status[1]
					|| d1.status[1] == d2.status[0]
				)
			)) {
				bind_points.Fill(comb1.points);
				bind_queue.push(comb2);
			}
		}
	}

// while (!bind_queue.empty()) {
// 	TrackBindingCombination comb = bind_queue.top();
// 	bind_queue.pop();
// 	std::cout << "========================\n"
// 		<< "points " << std::dec << comb.points << "\n"
// 		<< std::hex << "flag " << comb.flag[0] << " " << comb.flag[1] << "\n"
// 		<< "---------------first-----------------\n"
// 		<< "energy " << comb.energy[0][0] << " " << comb.energy[0][1] << "\n"
// 		<< "time " << comb.time[0][0] << " " << comb.time[0][1] << "\n"
// 		<< "x " << comb.x[0][0] << " " << comb.x[0][1] << "\n"
// 		<< "y " << comb.y[0][0] << " " << comb.y[0][1] << "\n"
// 		<< "status " << comb.status[0] << "\n"
// 		<< "---------------second-----------------\n"
// 		<< "energy " << comb.energy[1][0] << " " << comb.energy[1][1] << "\n"
// 		<< "time " << comb.time[1][0] << " " << comb.time[1][1] << "\n"
// 		<< "x " << comb.x[1][0] << " " << comb.x[1][1] << "\n"
// 		<< "y " << comb.y[1][0] << " " << comb.y[1][1] << "\n"
// 		<< "status " << comb.status[1] << "\n";
// }
	// binding events for statistics
	int binding_events = 0;
	// fill t0 events
	while (!bind_queue.empty() && d1_flag != 0 && d2_flag != 0 && t0.num < 7) {
		TrackBindingCombination comb = bind_queue.top();
		bind_queue.pop();
// 		// if (comb.points < 20) continue;
		unsigned int comb_d1_flag = comb.flag[0][0] | comb.flag[0][1];
		unsigned int comb_d2_flag = comb.flag[1][0] | comb.flag[1][1];
		if ((d1_flag & comb_d1_flag) != comb_d1_flag) continue;
		if ((d2_flag & comb_d2_flag) != comb_d2_flag) continue;
// std::cout << "\n" << t0.num << " d1 flag " << std::hex << comb.flag[0] << " / " << d1_flag << "\n"
// 	<< comb.flag[1] << " / " << d2_flag << "\n";
		d1_flag &= ~comb_d1_flag;
		d2_flag &= ~comb_d2_flag;
		for (size_t i = 0; i < 2; ++i) {
			t0.layer[t0.num+i] = 2;
			t0.flag[t0.num+i] = 0x3;
			t0.energy[t0.num+i][0] = comb.energy[i][0];
			t0.energy[t0.num+i][1] = comb.energy[i][1];
			t0.time[t0.num+i][0] = comb.time[i][0];
			t0.time[t0.num+i][1] = comb.time[i][1];
			t0.x[t0.num+i][0] = comb.x[i][0];
			t0.x[t0.num+i][1] = comb.x[i][1];
			t0.y[t0.num+i][0] = comb.y[i][0];
			t0.y[t0.num+i][1] = comb.y[i][1];
			t0.z[t0.num+i][0] = 100.0;
			t0.z[t0.num+i][1] = 111.76;
			t0.status[t0.num+i] = comb.status[i];
			t0.points[t0.num+i] = comb.points;
			t0.dssd_flag[t0.num+i][0] = comb.flag[0][i];
			t0.dssd_flag[t0.num+i][1] = comb.flag[1][i];
		}
		t0.num += 2;
		binding_events += 2;
	}

	// // particle types identified in T0D1 and T0D2
	// std::vector<int> d1d2_types;
	// // identify particles
	// for (int i = 0; i < t0.num; ++i) {
	// 	d1d2_types.push_back(
	// 		D1D2ParticleIdentify(t0.energy[i][0], t0.energy[i][1], d1d2_cuts)
	// 	);
	// }
	// particle types identified in T0D1 and T0D2 tail
	std::vector<int> d1d2_tail_types;
	// identify particle tails
	for (int i = 0; i < t0.num; ++i) {
		d1d2_tail_types.push_back(
			D1D2ParticleIdentify(t0.energy[i][0], t0.energy[i][1], d1d2_tails)
		);
	}

	// priority queue for T0D3 events
	std::priority_queue<D3Combination> d3_queue;

	// search for combinations and evaluate points
	for (const auto &d3 : d3_combinations) {
		if (d3.num == 4) continue;
		for (int i = 0; i < t0.num; ++i) {
			// if (d1d2_types[i] > 0) continue;
			if (d1d2_tail_types[i] == 0) continue;
			for (int j = 0; j < d3.num; ++j) {
				if (fabs(t0.x[i][1]-d3.x[j]) > 4) continue;
				if (fabs(t0.y[i][1]-d3.y[j]) > 4) continue;
				if (fabs(t0.time[i][1]-d3.time[j]) > comb_time_threshold) continue;
				d1d3dt.Fill(t0.time[i][0] - d3.time[j]);
				d1d3dx.Fill(t0.x[i][0] - d3.x[j]);
				d1d3dy.Fill(t0.y[i][0] - d3.y[j]);
				d2d3dt.Fill(t0.time[i][1] - d3.time[j]);
				d2d3dx.Fill(t0.x[i][1] - d3.x[j]);
				d2d3dy.Fill(t0.y[i][1] - d3.y[j]);

				int points = 0;

				if (fabs(t0.x[i][1]-d3.x[j]) < 1.0) points += 3;
				else if (fabs(t0.x[i][1]-d3.x[j]) < 2.0) points += 2;
				else if (fabs(t0.x[i][1]-d3.x[j]) < 3.0) points += 1;

				if (fabs(t0.y[i][1]-d3.y[j]) < 1.0) points += 3;
				else if (fabs(t0.y[i][1]-d3.y[j]) < 2.0) points += 2;
				else if (fabs(t0.y[i][1]-d3.y[j]) < 3.0) points += 1;

				if (fabs(t0.time[i][1] - d3.time[j]) < 20.0) points += 5;
				else if (fabs(t0.time[i][1] - d3.time[j]) < 80.0) points += 3;
				else if (fabs(t0.time[i][1] - d3.time[j]) < 150.0) points += 2;
				else if (fabs(t0.time[i][1] - d3.time[j]) < 200.0) points += 1;
				else points -= 1;

				D3Combination comb{
					d3.energy[j],
					d3.time[j],
					d3.x[j],
					d3.y[j],
					d3.flag[j],
					d3.status[j],
					points,
					i
				};
				d3_points.Fill(comb.points);
				d3_queue.push(comb);
			}
		}
	}

	// fill T0 events
	while (!d3_queue.empty() && d3_flag != 0 && t0.num < 8) {
		D3Combination comb = d3_queue.top();
		d3_queue.pop();
		// if (comb.points < 7) continue;
		if ((d3_flag & comb.flag) != comb.flag) continue;
		d3_flag &= ~comb.flag;
		if (t0.layer[comb.d1d2_index] == 3) continue;
		t0.layer[comb.d1d2_index] = 3;
		t0.flag[comb.d1d2_index] = 0x7;
		t0.energy[comb.d1d2_index][2] = comb.energy;
		t0.time[comb.d1d2_index][2] = comb.time;
		t0.x[comb.d1d2_index][2] = comb.x;
		t0.y[comb.d1d2_index][2] = comb.y;
		t0.z[comb.d1d2_index][2] = 123.52;
		t0.status[comb.d1d2_index] += comb.status*256;
		t0.points[comb.d1d2_index] += comb.points;
		t0.dssd_flag[comb.d1d2_index][2] = comb.flag;
	}

	return binding_events;
}


int T0::MergeAndTrack() {
	// T0D1 file name
	TString d1_file_name = TString::Format(
		"%s%st0d1-result-%s%04u.root",
		kGenerateDataPath,
		kNormalizeDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// input file
	TFile ipf(d1_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< d1_file_name << " failed.\n";
		return -1;
	}
	// T0D2 file name
	TString d2_file_name = TString::Format(
		"%s%st0d2-result-%s%04u.root",
		kGenerateDataPath,
		kNormalizeDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// add friend
	ipt->AddFriend("d2=tree", d2_file_name);
	// T0D3 file name
	TString d3_file_name = TString::Format(
		"%s%st0d3-result-%s%04u.root",
		kGenerateDataPath,
		kNormalizeDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// add friend
	ipt->AddFriend("d3=tree", d3_file_name);
	// input d1 event
	DssdFundamentalEvent d1_event;
	// input d2 event
	DssdFundamentalEvent d2_event;
	// input d3 event
	DssdFundamentalEvent d3_event;
	// setup input branches
	d1_event.SetupInput(ipt);
	d2_event.SetupInput(ipt, "d2.");
	d3_event.SetupInput(ipt, "d3.");

	// output file name
	TString output_file_name = TString::Format(
		"%s%st0-telescope-%s%04u.root",
		kGenerateDataPath,
		kTelescopeDir,
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// 1D histogram of T0D1 front adjacent strip dt
	TH1F d1_front_adj_dt("d1fadt", "#Deltat", 1000, -500, 500);
	// 1D histogram of T0D1 back adjacent strip dt
	TH1F d1_back_adj_dt("d1badt", "#Deltat", 1000, -500, 500);
	// 1D histogram of T0D2 front adjacent strip dt
	TH1F d2_front_adj_dt("d2fadt", "#Deltat", 1000, -500, 500);
	// 1D histogram of T0D2 back adjacent strip dt
	TH1F d2_back_adj_dt("d2badt", "#Deltat", 1000, -500, 500);
	// 1D histogram of T0D3 front adjacent strip dt
	TH1F d3_front_adj_dt("d3fadt", "#Deltat", 1000, -500, 500);
	// 1D histogram of T0D3 back adjacent strip dt
	TH1F d3_back_adj_dt("d3badt", "#Deltat", 1000, -500, 500);
	// 1D histogram of T0D1 side energy difference
	TH1F d1_side_de("d1sde", "#DeltaE", 1000, -2000, 2000);
	// 1D histogram of T0D1 side time difference
	TH1F d1_side_dt("d1sdt", "#Deltat", 1000, -500, 500);
	// 1D histogram of T0D2 side energy difference
	TH1F d2_side_de("d2sde", "#DeltaE", 1000, -2000, 2000);
	// 1D histogram of T0D2 side time difference
	TH1F d2_side_dt("d2sdt", "#Deltat", 1000, -500, 500);
	// 1D histogram of T0D3 side energy difference
	TH1F d3_side_de("d3sde", "#DeltaE", 1000, -2000, 2000);
	// 1D histogram of T0D3 side time difference
	TH1F d3_side_dt("d3sdt", "#Deltat", 1000, -500, 500);
	// 1D histogram of T0D1D2 dx
	TH1F d1d2_dx("d1d2dx", "#Deltax", 100, -10, 10);
	// 1D histogram of T0D1D2 dy
	TH1F d1d2_dy("d1d2dy", "#Deltay", 100, -5, 5);
	// 1D histogram of T0D1D2 dt
	TH1F d1d2_dt("d1d2dt", "#Deltat", 1000, -500, 500);
	// 1D histogram of D1D2 tracking points
	TH1F d1d2_points("d12p", "D1D2 points", 100, 0, 100);
	// 1D histogram of D1D2 binding points
	TH1F bind_points("bp", "D1D2 binding points", 100, 0, 100);
	// 1D histogram of T0D1D3 dx
	TH1F d1d3_dx("d1d3dx", "#Deltax", 100, -10, 10);
	// 1D histogram of T0D1D3 dy
	TH1F d1d3_dy("d1d3dy", "#Deltay", 100, -5, 5);
	// 1D histogram of T0D1D3 dt
	TH1F d1d3_dt("d1d3dt", "#Deltat", 1000, -500, 500);
	// 1D histogram of T0D2D3 dx
	TH1F d2d3_dx("d2d3dx", "#Deltax", 100, -10, 10);
	// 1D histogram of T0D2D3 dy
	TH1F d2d3_dy("d2d3dy", "#Deltay", 100, -5, 5);
	// 1D histogram of T0D2D3 dt
	TH1F d2d3_dt("d2d3dt", "#Deltat", 1000, -500, 500);
	// 1D histogram of D3 tracking points
	TH1F d3_points("d3p", "D3 points", 100, 0, 100);
	// output tree
	TTree opt("tree", "merge and track");
	// output event
	T0Event t0_event;
	// setup output branches
	t0_event.SetupOutput(&opt);


	// T0D1-D2 cuts
	std::vector<ParticleCut> d1d2_cuts;
	d1d2_cuts.push_back({2, 4, ReadCut("t0-d1d2-4He")});
	// d1d2_cuts.push_back({3, 6, ReadCut("t0-d1d2-6Li")});
	d1d2_cuts.push_back({3, 7, ReadCut("t0-d1d2-7Li")});
	// d1d2_cuts.push_back({4, 7, ReadCut("t0-d1d2-7Be")});
	// d1d2_cuts.push_back({4, 9, ReadCut("t0-d1d2-9Be")});
	d1d2_cuts.push_back({4, 10, ReadCut("t0-d1d2-10Be")});
	// d1d2_cuts.push_back({5, 10, ReadCut("t0-d1d2-10B")});
	// d1d2_cuts.push_back({5, 11, ReadCut("t0-d1d2-11B")});
	// d1d2_cuts.push_back({5, 12, ReadCut("t0-d1d2-12B")});
	// d1d2_cuts.push_back({5, 13, ReadCut("t0-d1d2-13B")});
	// d1d2_cuts.push_back({6, 12, ReadCut("t0-d1d2-12C")});
	// d1d2_cuts.push_back({6, 13, ReadCut("t0-d1d2-13C")});
	d1d2_cuts.push_back({6, 14, ReadCut("t0-d1d2-14C")});

	std::vector<ParticleCut> d1d2_tails;
	d1d2_tails.push_back({2, 0, ReadCut("t0-d1d2-tail-He")});
	d1d2_tails.push_back({3, 0, ReadCut("t0-d1d2-tail-Li")});
	d1d2_tails.push_back({4, 0, ReadCut("t0-d1d2-tail-Be")});

	// statistics
	long long particle_num[8];
	for (size_t i = 0; i < 8; ++i) particle_num[i] = 0;
	long long binding_events = 0;

	// total number of entries
	long long entries = ipt->GetEntries();
	// long long entries = 10;
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Merging and tracking T0   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
// std::cout << "--------------------------------\n"
// 	<< "entry  " << entry << "\n"
// 	<< "--------------------------------\n";
		// get event
		ipt->GetEntry(entry);
		std::vector<StripGroup> d1_front;
		std::vector<StripGroup> d1_back;
		// search for T0D1 strip combinations
		SearchStripCombinations(
			d1_event, d1_front_adj_dt, d1_back_adj_dt, d1_front, d1_back
		);
		// search for T0D1 side combinations
		std::vector<SideCombination> d1_comb = SearchSideCombinations(
			d1_front, d1_back, 0, d1_side_de, d1_side_dt
		);

		std::vector<StripGroup> d2_front;
		std::vector<StripGroup> d2_back;
		// search for T0D2 strip combinations
		SearchStripCombinations(
			d2_event, d2_front_adj_dt, d2_back_adj_dt, d2_front, d2_back
		);
		// search for T0D2 side combinations
		std::vector<SideCombination> d2_comb = SearchSideCombinations(
			d2_front, d2_back, 1, d2_side_de, d2_side_dt
		);

		std::vector<StripGroup> d3_front;
		std::vector<StripGroup> d3_back;
		// search for T0D2 strip combinations
		SearchStripCombinations(
			d3_event, d3_front_adj_dt, d3_back_adj_dt, d3_front, d3_back
		);
		// search for T0D2 side combinations
		std::vector<SideCombination> d3_comb = SearchSideCombinations(
			d3_front, d3_back, 2, d3_side_de, d3_side_dt
		);

		// tracking
		binding_events += SearchTrackCombinations(
			d1_comb, d2_comb, d3_comb,
			// d1d2_cuts,
			d1d2_tails,
			((1 << d1_event.front_hit) | (1 << (d1_event.back_hit+8))) - 0x101,
			((1 << d2_event.front_hit) | (1 << (d2_event.back_hit+8))) - 0x101,
			((1 << d3_event.front_hit) | (1 << (d3_event.back_hit+8))) - 0x101,
			d1d2_dt, d1d2_dx, d1d2_dy,
			d1d3_dt, d1d3_dx, d1d3_dy,
			d2d3_dt, d2d3_dx, d2d3_dy,
			d1d2_points, bind_points, d3_points,
			t0_event
		);
		// track T0D3

		if (t0_event.num <= 8 && t0_event.num > 0) {
			++particle_num[t0_event.num-1];
		}
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	long long total = 0;
	for (size_t i = 0; i < 8; ++i) {
		total += particle_num[i];
	}
	for (size_t i = 0; i < 8; ++i) {
		std::cout << "particle " << i+1 << "  events "
			<< particle_num[i] << " / " << total << "  "
			<< double(particle_num[i]) / double(total) << "\n";
	}
	std::cout << "Add binding events " << binding_events << "\n";

	// save histograms
	d1_front_adj_dt.Write();
	d1_back_adj_dt.Write();
	d2_front_adj_dt.Write();
	d2_back_adj_dt.Write();
	d3_front_adj_dt.Write();
	d3_back_adj_dt.Write();
	d1_side_de.Write();
	d1_side_dt.Write();
	d2_side_de.Write();
	d2_side_dt.Write();
	d3_side_de.Write();
	d3_side_dt.Write();
	d1d2_dx.Write();
	d1d2_dy.Write();
	d1d2_dt.Write();
	d1d3_dx.Write();
	d1d3_dy.Write();
	d1d3_dt.Write();
	d2d3_dx.Write();
	d2d3_dy.Write();
	d2d3_dt.Write();
	d1d2_points.Write();
	bind_points.Write();
	d3_points.Write();
	// save tree
	opt.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}


}		// namespace ribll