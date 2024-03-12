#include <iostream>

#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TRandom3.h>
#include <TString.h>

#include "include/event/generate_event.h"
#include "include/event/dssd_event.h"
#include "include/event/ssd_event.h"
#include "include/calculator/range_energy_calculator.h"

using namespace ribll;

// 能量分辨（sigma）模拟目标（道址）
// T0D1-120(be),130(he)
// T0D2-210(be),150(he)
// T0D3-30(be),45(he)
constexpr double t0_energy_sigma[7] = {0.5, 0.8, 0.16, 0.1, 0.1, 0.1, 1.0};
// T0 DSSD center position for correcting
constexpr double t0_center[3][2] = {
	{0.0, 0.0},
	{-1.03, -0.86},
	{-0.95, -0.8}
};
// T0 detect threshold
constexpr double t0_energy_threshold[6] = {
	1000.0, 1000.0, 1000.0, 2000.0, 500.0, 800.0
};
// T0 corrected thickness
constexpr double thickness[6] = {
	1010.0, 1504.0, 1501.0, 1480.0, 1500.0, 1650.0
};

int main(int argc, char **argv) {
	int run = 0;
	if (argc > 1) {
		run = atoi(argv[1]);
	}
	if (run < 0 || run > 1) {
		std::cout << "Usage: " << argv[0] << "[run]\n"
			<< "  run        run number, default is 0\n";
	}

	// input generate data file name
	TString generate_file_name = TString::Format(
		"%s%sgenerate-%04d.root",
		kGenerateDataPath,
		kSimulateDir,
		run
	);
	// generate data file
	TFile generate_file(generate_file_name, "read");
	// input generate tree
	TTree *generate_tree = (TTree*)generate_file.Get("tree");
	if (!generate_tree) {
		std::cerr << "Error: Get tree from "
			<< generate_file_name << " failed.\n";
		return -1;
	}
	// generate event
	GenerateEvent event;
	// setup input branches
	event.SetupInput(generate_tree);

	// output DSSD files name
	TString dssd_file_names[3];
	for (int i = 0; i < 3; ++i) {
	 	dssd_file_names[i].Form(
			"%s%st0d%d-result-sim-ta-%04d.root",
			kGenerateDataPath,
			kNormalizeDir,
			i+1,
			run
		);
	}
	// output DSSD file
	TFile *dssd_files[3];
	for (int i = 0; i < 3; ++i) {
		dssd_files[i] = new TFile(dssd_file_names[i], "recreate");
	}
	// output DSSD trees
	TTree *dssd_trees[3];
	for (int i = 0; i < 3; ++i) {
		dssd_files[i]->cd();
		dssd_trees[i] = new TTree("tree", "simulated DSSD data");
	}
	// output DSSD events
	DssdFundamentalEvent dssd_events[3];
	// setup output branches
	for (int i = 0; i < 3; ++i) {
		dssd_events[i].SetupOutput(dssd_trees[i]);
	}

	// output SSD file names
	TString ssd_file_names[3];
	for (int i = 0; i < 3; ++i) {
		ssd_file_names[i].Form(
			"%s%st0s%d-fundamental-sim-ta-%04d.root",
			kGenerateDataPath,
			kFundamentalDir,
			i+1,
			run
		);
	}
	// output SSD file
	TFile *ssd_files[3];
	for (int i = 0; i < 3; ++i) {
		ssd_files[i] = new TFile(ssd_file_names[i], "recreate");
	}
	// output SSD tree
	TTree *ssd_trees[3];
	for (int i = 0; i < 3; ++i) {
		ssd_files[i]->cd();
		ssd_trees[i] = new TTree("tree", "simulated SSD data");
	}
	// output SSD event
	SsdEvent ssd_events[3];
	// setup output branches
	for (int i = 0; i < 3; ++i) {
		ssd_events[i].SetupOutput(ssd_trees[i]);
	}

	// output detect file for monitoring
	// output detect file name
	TString detect_file_name = TString::Format(
		"%s%st0-detect-%04d.root",
		kGenerateDataPath, kSimulateDir, run
	);
	// output file
	TFile detect_file(detect_file_name, "recreate");
	// histogram of comparing T0Dx front-back energy
	TH1F *hist_fb_energy[3];
	for (int i = 0; i < 3; ++i) {
		hist_fb_energy[i] = new TH1F(
			TString::Format("hd%dde", i+1),
			TString::Format("T0D%d front back energy", i+1),
			1000, -1000, 1000
		);
	}
	// T0 detect tree
	TTree opt("tree", "T0 detect information");
	// output position DSSD flag
	// -4: 1 particle, cut by energy threshold
	// -3: 2 particle, cut by energy threshold
	// -2: out of range
	// -1: cannot reach this layer
	// 0: nothing, invalid
	// 1: 1 particle
	// 2: 2 seperated strip in both sides
	// 3: 2 seperated strip in front side, 2 adjacent strips in back side
	// 4: 2 seperated strip in front side, 1 same strip in back side
	// 5: 2 adjacent strips in front side, 2 seperated strip in back side
	// 6: 2 adjacent strips in both sides
	// 7: 2 adjacent strips in front side, 1 same strip in back side
	// 8: 1 same strip in front side, 2 seperated strip in back side
	// 9: 1 same strip in front side, 2 adjacent strips in back side
	// 10: 1 same strip in both sides
	int dssd_flag[3];

	// output position T0 flag
	// 0: impossible t0 rebuild
	// 1: only 2 strips events
	// 2: 2 strips events, includes adjacent strips
	// 3: includes 1 strip in one side, 2 seperated strip in the other side
	// 4: includes 1 strip in one side, 2 adjacent strips in the other side
	// 5: includes 1 strip in both side
	int t0_flag;

	// times cut by energy threshold
	int cut_time;
	// stopped layer: 0-T0D1, 1-T0D2, 2-T0D3, 3-T0S1, 4-T0S2, 5-T0S3, 6-T0CsI
	int t0_layer[2];
	double t0_lost_energy[2][7];
	double t0_front_energy[2][3];
	double t0_back_energy[2][3];
	double t0_front_channel[2][3];
	double t0_back_channel[2][3];
	double t0_ssd_energy[2][3];
	double t0_ssd_channel[2][3];

	// setup output branches
	opt.Branch("dssd_flag", dssd_flag, "df[3]/I");
	opt.Branch("t0_flag", &t0_flag, "tf/I");
	opt.Branch("cut_time", &cut_time, "ct/I");
	opt.Branch("layer", t0_layer, "layer[2]/I");
	opt.Branch("lost_energy", t0_lost_energy, "le[2][7]/D");
	opt.Branch("front_energy", t0_front_energy, "fe[2][3]/D");
	opt.Branch("back_energy", t0_back_energy, "be[2][3]/D");
	opt.Branch("front_channel", t0_front_channel, "fc[2][3]/D");
	opt.Branch("back_channel", t0_back_channel, "bc[2][3]/D");
	opt.Branch("ssd_energy", t0_ssd_energy, "se[2][3]/D");
	opt.Branch("ssd_channel", t0_ssd_channel, "sc[2][3]/D");


	// initialize random number generator
	TRandom3 generator(0);

	// initialize calculators
	elc::RangeEnergyCalculator be10_calculator("10Be", "Si");
	elc::RangeEnergyCalculator he4_calculator("4He", "Si");
	elc::RangeEnergyCalculator* frag_calculators[2] = {
		&be10_calculator, &he4_calculator
	};

	// statistics
	// total counts
	int total_count[3] = {0, 0, 0};
	// out of range counts
	int out_of_range_count[3] = {0, 0, 0};
	// cut by energy threshold counts
	int p1_under_threshold_count[3] = {0, 0, 0};
	int p2_under_threshold_count[3] = {0, 0, 0};
	// invalid counts, should be zero
	int invalid_count[3] = {0, 0, 0};
	// different kinds strip distribution counts
	int valid_counts[3][10];
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 10; ++j) valid_counts[i][j] = 0;
	}
	// total valid counts for T0
	int t0_total_count = 0;
	// only 2 strip events in T0
	int t0_only_two_count = 0;
	// 2 strip events includes adjacent strips in both sides
	int t0_both_adjdacent_count = 0;
	// includes 1-2 events, flag 4 and 8 events
	// 1 same strip in one side, 2 seperated strips in the other side
	int t0_one_two_count = 0;
	// includes 1-2a events, flag 7 and 9 events
	// 1 same strip in one side, 2 adjacent strip2 in the other side
	int t0_one_adjacent_count = 0;
	// includes 1-1 events, flag 10 events
	// 1 strip in both sides
	int t0_both_one_count = 0;



	// initialize useless variables
	for (int i = 0; i < 3; ++i) {
		dssd_events[i].cfd_flag = 0;
		for (int j = 0; j < 8; ++j) {
			dssd_events[i].front_time[j] = 0;
			dssd_events[i].back_time[j] = 0;
			dssd_events[i].front_decode_entry[j] = 0;
			dssd_events[i].back_decode_entry[j] = 0;
			dssd_events[i].front_fundamental_index[j] = 0;
			dssd_events[i].back_fundamental_index[j] = 0;
		}
		ssd_events[i].cfd_flag = 0;
	}

	// show start
	printf("Simulating T0 detect   0%%");
	// total entries, fow showing process
	long long entries = generate_tree->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100;
	// detecting simulation
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get entry
		generate_tree->GetEntry(entry);

		// initialize
		for (int i = 0; i < 3; ++i) {
			dssd_events[i].front_hit = 0;
			dssd_events[i].back_hit = 0;
			ssd_events[i].time = -1e5;
			ssd_events[i].energy = 0.0;
		}
		for (int i = 0; i < 3; ++i) dssd_flag[i] = 0;
		t0_flag = 0;
		cut_time = 0;

		// calculate t0 energy
		double range[2];
		double residual_energy[2];
		for (size_t i = 0; i < 2; ++i) {
			for (size_t j = 0; j < 7; ++j) {
				t0_lost_energy[i][j] = 0.0;
			}
			range[i] = frag_calculators[i]->Range(event.fragment_kinetic[i]);
			residual_energy[i] = event.fragment_kinetic[i];

			for (size_t j = 0; j < 6; ++j) {
				double thick = thickness[j] / cos(event.fragment_theta[i]);
				if (range[i] < thick) {
					t0_layer[i] = j;
					t0_lost_energy[i][j] = residual_energy[i];
					break;
				} else {
					t0_lost_energy[i][j] =
						residual_energy[i] - frag_calculators[i]->Energy(range[i] - thick);
					range[i] -= thick;
					residual_energy[i] -= t0_lost_energy[i][j];
					if (j == 5) {
						t0_layer[i] = 6;
						t0_lost_energy[i][6] = residual_energy[i];
					}
				}
			}
		}

		// check layers
		int stop_t0d2 = 0;
		if (t0_layer[0] > 0) ++stop_t0d2;
		if (t0_layer[1] > 0) ++stop_t0d2;
		if (stop_t0d2 == 0) dssd_flag[1] = -1;
		else if (stop_t0d2 == 1) dssd_flag[1] = 1;

		int stop_t0d3 = 0;
		if (t0_layer[0] > 1) ++stop_t0d3;
		if (t0_layer[1] > 1) ++stop_t0d3;
		if (stop_t0d3 == 0) dssd_flag[2] = -1;
		else if (stop_t0d3 == 1) dssd_flag[2] = 1;

		// simulate detected energy
		for (int i = 0; i < 2; ++i) {
			// DSSD energy
			for (int j = 0; j < 3; ++j) {
				t0_front_energy[i][j] =
					t0_lost_energy[i][j]
					+ generator.Gaus(0.0, t0_energy_sigma[j]);
				t0_back_energy[i][j] =
					t0_lost_energy[i][j]
					+ generator.Gaus(0.0, t0_energy_sigma[j]);
			}
			// SSD energy
			for (int j = 0; j < 3; ++j) {
				if (t0_lost_energy[i][j+3] == 0.0) continue;
				t0_ssd_energy[i][j] =
					t0_lost_energy[i][j+3]
					+ generator.Gaus(0.0, t0_energy_sigma[j+3]);
			}
		}

		// convert energy to channel
		for (int i = 0; i < 2; ++i) {
			// DSSD channel
			for (int j = 0; j < 3; ++j) {
				t0_front_channel[i][j] =
					(t0_front_energy[i][j] - t0_param[j][0])
					/ t0_param[j][1];
				t0_back_channel[i][j] =
					(t0_back_energy[i][j] - t0_param[j][0])
					/ t0_param[j][1];
			}
			// SSD channel
			for (int j = 0; j < 3; ++j) {
				t0_ssd_channel[i][j] =
					(t0_ssd_energy[i][j] - t0_param[j+3][0])
					/ t0_param[j+3][1];
			}
		}

		// calculate T0 position
		for (size_t i = 0; i < 2; ++i) {
			// check range
			if (event.fragment_x[i] > -32.0 && event.fragment_x[i] < 32.0) {
				// T0D1 position
				int t0d1_x_strip = int(event.fragment_x[i] + 32.0);
				int t0d1_y_strip = int(event.fragment_y[i] + 32.0);

				// fill T0D1 front information
				bool t0d1_front_fill = false;
				for (int j = 0; j < dssd_events[0].front_hit; ++j) {
					if (dssd_events[0].front_strip[j] == t0d1_y_strip) {
						dssd_events[0].front_energy[j] += t0_front_channel[i][0];
						t0d1_front_fill = true;
						break;
					}
				}
				if (!t0d1_front_fill) {
					int &j = dssd_events[0].front_hit;
					dssd_events[0].front_strip[j] = t0d1_y_strip;
					dssd_events[0].front_energy[j] = t0_front_channel[i][0];
					++j;
				}
				// fill T0D1 back information
				bool t0d1_back_fill = false;
				for (int j = 0; j < dssd_events[0].back_hit; ++j) {
					if (dssd_events[0].back_strip[j] == t0d1_x_strip) {
						dssd_events[0].back_energy[j] += t0_back_channel[i][0];
						t0d1_back_fill = true;
						break;
					}
				}
				if (!t0d1_back_fill) {
					int &j = dssd_events[0].back_hit;
					dssd_events[0].back_strip[j] = t0d1_x_strip;
					dssd_events[0].back_energy[j] = t0_back_channel[i][0];
					++j;
				}
			} else {
				dssd_flag[0] = -2;
			}


			// T0D2, D3 position
			for (int j = 1; j < 3; ++j) {
				// check layer
				if (t0_layer[i] < j) break;
				// get position
				double z = j == 1 ? 111.76 : 123.52;
				double r = z * tan(event.fragment_theta[i]);
				double x = r * cos(event.fragment_phi[i]) + event.target_x;
				double y = r * sin(event.fragment_phi[i]) + event.target_y;
				// check range
				if (x <= -32.0 || x >= 32.0 || y <= -32.0 || y >= 32.0) {
					dssd_flag[j] = -2;
					continue;
				}
				// position correct
				x -= t0_center[j][0];
				y -= t0_center[j][1];
				// strip number
				int xstrip = int((x+32.0)/2.0);
				int ystrip = int((y+32.0)/2.0);

				// fill T0Dx front information
				bool front_fill = false;
				for (int k = 0; k < dssd_events[j].front_hit; ++k) {
					if (dssd_events[j].front_strip[k] == xstrip) {
						dssd_events[j].front_energy[k] += t0_front_channel[i][j];
						front_fill = true;
						break;
					}
				}
				if (!front_fill) {
					int &k = dssd_events[j].front_hit;
					dssd_events[j].front_strip[k] = xstrip;
					dssd_events[j].front_energy[k] = t0_front_channel[i][j];
					++k;
				}

				// fill T0Dx back information
				bool back_fill = false;
				for (int k = 0; k < dssd_events[j].back_hit; ++k) {
					if (dssd_events[j].back_strip[k] == ystrip) {
						dssd_events[j].back_energy[k] += t0_back_channel[i][j];
						back_fill = true;
						break;
					}
				}
				if (!back_fill) {
					int &k = dssd_events[j].back_hit;
					dssd_events[j].back_strip[k] = ystrip;
					dssd_events[j].back_energy[k] = t0_back_channel[i][j];
					++k;
				}
			}
		}

		// fill SSD information
		for (int i = 0; i < 2; ++i) {
			for (int j = 3; j < 6; ++j) {
				if (t0_layer[i] < j) break;
				ssd_events[j-3].time = 0.0;
				ssd_events[j-3].energy += t0_ssd_channel[i][j-3];
			}
		}


		// sort DSSD events
		for (int i = 0; i < 3; ++i) {
			dssd_events[i].Sort();
		}


		// check threshold
		// check DSSD
		for (int i = 0; i < 3; ++i) {
			// check front side
			for (int j = dssd_events[i].front_hit-1; j >= 0; --j) {
				if (dssd_events[i].front_energy[j] < t0_energy_threshold[i]) {
					--dssd_events[i].front_hit;
					if (dssd_flag[i] == 0) dssd_flag[i] = -3;
					else if (dssd_flag[i] == 1) dssd_flag[i] = -4;
					++cut_time;
				}
			}
			for (int j = dssd_events[i].back_hit-1; j >= 0; --j) {
				if (dssd_events[i].back_energy[j] < t0_energy_threshold[i]) {
					--dssd_events[i].back_hit;
					if (dssd_flag[i] == 0) dssd_flag[i] = -3;
					else if (dssd_flag[i] == 1) dssd_flag[i] = -4;
				}
			}
		}


		// check SSD
		for (int i = 0; i < 3; ++i) {
			if (
				ssd_events[i].time > -9e4
				&& ssd_events[i].energy < t0_energy_threshold[i+3]
			) {
				ssd_events[i].time = -1e5;
				ssd_events[i].energy = 0.0;
				++cut_time;
			}
		}

// std::cout << "layer " << t0_layer[0] << ", " << t0_layer[1] << "\n"
// 	<< "T0D1 hit " << dssd_events[0].front_hit << ", " << dssd_events[0].back_hit << "\n"
// 	<< "T0D2 hit " << dssd_events[1].front_hit << ", " << dssd_events[1].back_hit << "\n"
// 	<< "T0D3 hit " << dssd_events[2].front_hit << ", " << dssd_events[2].back_hit << "\n"
// 	<< "Flag " << dssd_flag[0] << ", " << dssd_flag[1] << ", " << dssd_flag[2] << "\n"
// 	<< "Stop number " << stop_t0d2 << ", " << stop_t0d3 << "\n";

		// get flag
		for (int i = 0; i < 3; ++i) {
			if (dssd_flag[i] != 0) continue;
			// front side flag
			int front_flag = -1;
			// check front side
			if (dssd_events[i].front_hit == 1) front_flag = 2;
			else if (dssd_events[i].front_hit == 2) {
				// check if it's adjacent strip
				int delta_strip = abs(
					dssd_events[i].front_strip[0]
					- dssd_events[i].front_strip[1]
				);

				if (delta_strip == 1) front_flag = 1;
				else front_flag = 0;
			} else {
				std::cerr << "Error: Contradict T0 flag and front hit: entry "
					<< entry << " T0D" << i+1 << "\n";
				return -1;
			}
			// back side flag
			int back_flag = -1;
			// check back side flag
			if (dssd_events[i].back_hit == 1) back_flag = 2;
			else if (dssd_events[i].back_hit == 2) {
				// check if it's adjacent strip
				int delta_strip = abs(
					dssd_events[i].back_strip[0]
					- dssd_events[i].back_strip[1]
				);

				if (delta_strip == 1) back_flag = 1;
				else back_flag = 0;
			} else {
				std::cerr << "Error: Contradict T0 flag and back hit: entry "
					<< entry << " T0D" << i+1 << "\n";
				return -1;
			}
			dssd_flag[i] = front_flag*3 + back_flag + 2;
		}

		// count flags
		for (int i = 0; i < 3; ++i) {
			if (dssd_flag[i] != 0 && dssd_flag[i] != -1) ++total_count[i];
			if (dssd_flag[i] == -4) ++p1_under_threshold_count[i];
			else if (dssd_flag[i] == -3) ++p2_under_threshold_count[i];
			else if (dssd_flag[i] == -2) ++out_of_range_count[i];
			else if (dssd_flag[i] == 0) ++invalid_count[i];
			else if (dssd_flag[i] >= 1 && dssd_flag[i] <= 10) {
				++valid_counts[i][dssd_flag[i]-1];
			}
		}
		// count T0 statistics
		if (
			dssd_flag[0] >= 2 && dssd_flag[0] <= 10
			&& dssd_flag[1] >= 2 && dssd_flag[1] <= 10
			&& (
				(dssd_flag[2] >= 1 && dssd_flag[2] <= 10)
				|| dssd_flag[2] == -1
			)
		) {
			++t0_total_count;
			// shorter dssd flag
			int short_flag[3] = {0, 0, 0};
			for (int i = 0; i < 3; ++i) {
				if (dssd_flag[i] == 1) {
					short_flag[i] = 1;
				} else if (dssd_flag[i] == 5) {
					short_flag[i] = 3;
				} else if (dssd_flag[i] == 4 || dssd_flag[i] == 8) {
					short_flag[i] = 4;
				} else if (dssd_flag[i] == 7 || dssd_flag[i] == 9) {
					short_flag[i] = 5;
				} else if (dssd_flag[i] == 10) {
					short_flag[i] = 6;
				} else {
					short_flag[i] = 2;
				}
			}
			t0_flag = -1;
			// get T0 flag
			if (
				short_flag[0] == 2
				&& short_flag[1] == 2
				&& short_flag[2] <= 2
			) {
				t0_flag = 1;
				++t0_only_two_count;
			} else if (
				(	short_flag[0] == 3
					|| short_flag[1] == 3
					|| short_flag[2] == 3
				) && (
					short_flag[0] <= 3
					&& short_flag[1] <= 3
					&& short_flag[2] <= 3
				)
			) {
				t0_flag = 2;
				++t0_both_adjdacent_count;
			} else if (
				(
					short_flag[0] == 4
					|| short_flag[1] == 4
					|| short_flag[2] == 4
				) && (
					short_flag[0] <= 4
					&& short_flag[1] <= 4
					&& short_flag[2] <= 4
				)
			) {
				t0_flag = 3;
				++t0_one_two_count;
			} else if (
				(
					short_flag[0] == 5
					|| short_flag[1] == 5
					|| short_flag[2] == 5
				) && (
					short_flag[0] <= 5
					&& short_flag[1] <= 5
					&& short_flag[2] <= 5
				)
			) {
				t0_flag = 4;
				++t0_one_adjacent_count;
			} else if (
				short_flag[0] == 6
				|| short_flag[1] == 6
				|| short_flag[2] == 6
			) {
				t0_flag = 5;
				++t0_both_one_count;
			} else {
				std::cerr << "Error: Code error, invalid short_flag.\n";
				return -1;
			}

			if (t0_flag == -1) {
				std::cerr << "Error: Code error, invalid T0 flag.\n";
				return -1;
			}

		} else {
			t0_flag = 0;
		}

		// fill to histogram
		for (int i = 0; i < 3; ++i) {
			if (
				dssd_events[i].front_hit == 2
				&& dssd_events[i].back_hit == 2
			) {
				hist_fb_energy[i]->Fill(
					dssd_events[i].front_energy[0]
					- dssd_events[i].back_energy[0]
				);
				hist_fb_energy[i]->Fill(
					dssd_events[i].front_energy[1]
					- dssd_events[i].back_energy[1]
				);
			}
		}

		// fill DSSD tree
		for (int i = 0; i < 3; ++i) {
			dssd_trees[i]->Fill();
		}
		// fill SSD tree
		for (int i = 0; i < 3; ++i) {
			ssd_trees[i]->Fill();
		}
		// fill T0 detect output tree
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// fit histograms
	for (int i = 0; i < 3; ++i) {
		TF1 *f1 = new TF1(TString::Format("f%d", i+1), "gaus", -500, 500);
		f1->SetParameter(1, 0.0);
		f1->SetParameter(2, 100.0);
		hist_fb_energy[i]->Fit(f1, "RQ+");
		std::cout << "T0D" << i+1 << ": " << f1->GetParameter(2) << "\n";
	}

	// save histograms and tree
	detect_file.cd();
	for (int i = 0; i < 3; ++i) {
		hist_fb_energy[i]->Write();
	}
	opt.Write();
	detect_file.Close();

	// save trees and close output files
	for (int i = 0; i < 3; ++i) {
		dssd_files[i]->cd();
		dssd_trees[i]->Write();
		dssd_files[i]->Close();
	}
	for (int i = 0; i < 3; ++i) {
		ssd_files[i]->cd();
		ssd_trees[i]->Write();
		ssd_files[i]->Close();
	}
	// close input file
	generate_file.Close();

	// output T0D2 pixel file name
	TString t0d2_pixel_file_name = TString::Format(
		"%s%sshow-t0d2-pixel-%04d.root", kGenerateDataPath, kShowDir, run
	);
	// output T0D2 pixel file
	TFile t0d2_pixel_file(t0d2_pixel_file_name, "recreate");
	// 2D graph of beam energy resolution
	TH2F hist_beam_resolution(
		"hbr", "beam energy resolution", 32, 0, 32, 32, 0, 32
	);
	// save T0D2 resolution
	t0d2_pixel_file.cd();
	hist_beam_resolution.Write();
	t0d2_pixel_file.Close();

	// show statistics
	for (int i = 0; i < 3; ++i) {
		double total = double(total_count[i]);
		std::cout << "---------- T0D" << i+1 << " ----------\n"
			<< "Invalid " << invalid_count[i] << "\n"
			<< "Total " << total_count[i] << "\n"
			<< "1 particle under threshold " << p1_under_threshold_count[i]
			<< ", " << p1_under_threshold_count[i] / total << "\n"
			<< "2 particle under threshold " << p2_under_threshold_count[i]
			<< ", " << p2_under_threshold_count[i] / total << "\n"
			<< "Out of range " << out_of_range_count[i]
			<< ", " << out_of_range_count[i] / total << "\n"
			<< "1 particle " << valid_counts[i][0]
			<< ", " << valid_counts[i][0] / total << "\n"
			<< "Front 2 back 2 " << valid_counts[i][1]
			<< ", " << valid_counts[i][1] / total << "\n"
			<< "Front 2 back 2a " << valid_counts[i][2]
			<< ", " << valid_counts[i][2] / total << "\n"
			<< "Front 2 back 1 " << valid_counts[i][3]
			<< ", " << valid_counts[i][3] / total << "\n"
			<< "Front 2a back 2 " << valid_counts[i][4]
			<< ", " << valid_counts[i][4] / total << "\n"
			<< "Front 2a back 2a " << valid_counts[i][5]
			<< ", " << valid_counts[i][5] / total << "\n"
			<< "Front 2a back 1 " << valid_counts[i][6]
			<< ", " << valid_counts[i][6] / total << "\n"
			<< "Front 1 back 2 " << valid_counts[i][7]
			<< ", " << valid_counts[i][7] / total << "\n"
			<< "Front 1 back 2a " << valid_counts[i][8]
			<< ", " << valid_counts[i][8] / total << "\n"
			<< "Front 1 back 1 " << valid_counts[i][9]
			<< ", " << valid_counts[i][9] / total << "\n";
	}
	double t0_total = double(t0_total_count);
	std::cout << "---------- T0" << " ----------\n"
		<< "Total " << t0_total_count << "\n"
		<< "Only 2 " << t0_only_two_count << ", "
		<< t0_only_two_count / t0_total << "\n"
		<< "Include adjacent-adjacent " << t0_both_adjdacent_count << ", "
		<< t0_both_adjdacent_count / t0_total << "\n"
		<< "Include 1-2 " << t0_one_two_count << ", "
		<< t0_one_two_count / t0_total << "\n"
		<< "Include 1-2a " << t0_one_adjacent_count << ", "
		<< t0_one_adjacent_count / t0_total << "\n"
		<< "Include 1-1 " << t0_both_one_count << ", "
		<< t0_both_one_count / t0_total << "\n";

	return 0;
}