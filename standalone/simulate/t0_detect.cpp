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


int main() {
	// input generate data file name
	TString generate_file_name = TString::Format(
		"%s%sgenerate.root",
		kGenerateDataPath,
		kSimulateDir
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
			"%s%st0d%d-result-sim-ta-0000.root",
			kGenerateDataPath,
			kNormalizeDir,
			i+1
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
			"%s%st0s%d-fundamental-sim-ta-0000.root",
			kGenerateDataPath,
			kFundamentalDir,
			i+1
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
		"%s%st0-detect.root", kGenerateDataPath, kSimulateDir
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

	// initialize random number generator
	TRandom3 generator(0);

	// initialize calculators
	elc::RangeEnergyCalculator be10_calculator("10Be", "Si");
	elc::RangeEnergyCalculator he4_calculator("4He", "Si");
	elc::RangeEnergyCalculator* frag_calculators[2] = {
		&be10_calculator, &he4_calculator
	};

	// temporary variables
	double t0_lost_energy[2][7];
	double t0_front_energy[2][3];
	double t0_back_energy[2][3];
	double t0_front_channel[2][3];
	double t0_back_channel[2][3];
	double t0_ssd_energy[2][3];
	double t0_ssd_channel[2][3];
	double t0_layer[2];

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
				double thick = t0_thickness[j] / cos(event.fragment_theta[i]);
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
			}


			// T0D2, D3 position
			for (size_t j = 1; j < 3; ++j) {
				// check layer
				if (t0_layer[i] < j) break;
				// get position
				double z = j == 1 ? 111.76 : 123.52;
				double r = z * tan(event.fragment_theta[i]);
				double x = r * cos(event.fragment_phi[i]) + event.target_x;
				double y = r * sin(event.fragment_phi[i]) + event.target_y;
				// check range
				if (x <= -32.0 || x >= 32.0) continue;
				if (y <= -32.0 || y >= 32.0) continue;
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

		// check threshold
		// check DSSD
		for (int i = 0; i < 3; ++i) {
			// check front side
			bool front_change = true;
			while (front_change) {
				front_change = false;
				// for convenient
				int &fhit = dssd_events[i].front_hit;
				unsigned short *fs = dssd_events[i].front_strip;
				double *fe = dssd_events[i].front_energy;
				for (int j = 0; j < fhit; ++j) {
					if (fe[j] < t0_energy_threshold[i]) {
						fs[j] = fs[fhit-1];
						fe[j] = fe[fhit-1];
						--fhit;
						front_change = true;
						break;
					}
				}
			}
			// check back side
			bool back_change = true;
			while (back_change) {
				back_change = false;
				// for convenient
				int &bhit = dssd_events[i].back_hit;
				unsigned short *bs = dssd_events[i].back_strip;
				double *be = dssd_events[i].back_energy;
				for (int j = 0; j < bhit; ++j) {
					if (be[j] < t0_energy_threshold[i]) {
						bs[j] = bs[bhit-1];
						be[j] = be[bhit-1];
						--bhit;
						back_change = true;
						break;
					}
				}
			}
		}
		// check SSD
		for (int i = 0; i < 3; ++i) {
			if (
				ssd_events[i].time < -9e4
				&& ssd_events[i].energy < t0_energy_threshold[i+3]
			) {
				ssd_events[i].time = -1e5;
				ssd_events[i].energy = 0.0;
			}
		}


		// sort events
		for (int i = 0; i < 3; ++i) {
			dssd_events[i].Sort();
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

	// save histograms
	detect_file.cd();
	for (int i = 0; i < 3; ++i) {
		hist_fb_energy[i]->Write();
	}
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
		"%s%sshow-t0d2-pixel-0000.root", kGenerateDataPath, kShowDir
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

	return 0;
}