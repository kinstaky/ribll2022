#include "include/telescope/t0.h"

#include <Math/Vector3D.h>
#include <TH1F.h>

#include "include/event/dssd_event.h"
#include "include/event/t0_event.h"
#include "include/statistics/track_statistics.h"
#include "include/statistics/track_dssd_statistics.h"

namespace ribll {

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
	// d1 event
	DssdMergeEvent d1_event;
	d1_event.SetupInput(ipt);
	// d2 event
	DssdMergeEvent d2_event;
	d2_event.SetupInput(ipt, "d2.");
	// d3 event
	DssdMergeEvent d3_event;
	d3_event.SetupInput(ipt, "d3.");

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

}		// namespace ribll