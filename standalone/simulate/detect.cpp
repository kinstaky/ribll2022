#include <cmath>
#include <iostream>

#include <TFile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/generate_event.h"
#include "include/event/detect_event.h"
#include "include/calculator/range_energy_calculator.h"

using namespace ribll;

// parameters control the affect
constexpr bool t0_strips_affect = true;
constexpr bool taf_strips_affect = true;
constexpr bool ppac_strips_affect = true;
constexpr bool t0_energy_affect = true;
constexpr bool taf_energy_affect = true;

constexpr double pi = 3.1415926;
constexpr double u = 931.494;
constexpr double c14_mass = 13.9999505089 * u;
constexpr double be10_mass = 10.0113403769 * u;
constexpr double he4_mass = 4.0015060943 * u;
constexpr double h2_mass = 2.0135531980 * u;

constexpr double ppac_xz[3] = {-695.2, -454.2, -275.2};
constexpr double ppac_yz[3] = {-689.2, -448.2, -269.2};

constexpr double state_start[3] = {12.5, 16.0, 18.5};

double SimpleFit(const double *x, double *y, double &k, double &b) {
	int n = 3;
	double sumx = 0.0;
	double sumy = 0.0;
	double sumxy = 0.0;
	double sumx2 = 0.0;
	for (int i = 0; i < n; ++i) {
		sumx += x[i];
		sumy += y[i];
		sumxy += x[i] * y[i];
		sumx2 += x[i] * x[i];
	}
	k = (sumxy - sumx*sumy/double(n)) / (sumx2 - sumx*sumx/double(n));
	b = (sumy - k*sumx) / double(n);
	double chi2 = 0.0;
	for (int i = 0; i < n; ++i) {
		double t = y[i] - k*x[i] - b;
		chi2 += t * t;
	}
	return chi2;
}

int main() {
	// generate data file name
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

	// output file name
	TString detect_file_name = TString::Format(
		"%s%sdetect.root",
		kGenerateDataPath,
		kSimulateDir
	);
	// output file
	TFile detect_file(detect_file_name, "recreate");
	// fake q value histogram
	TH1F hist_fake_q("hfq", "fake q value", 1000, -20, -10);
	// T0D1,D2 x offset
	TH1F hist_d1d2_x_offset("hd1d2xo", "D1D2 x offset", 20, -10, 10);
	// T0D1,D2 x offset
	TH1F hist_d1d2_y_offset("hd1d2yo", "D1D2 y offset", 20, -10, 10);
	// T0D1,D3 x offset
	TH1F hist_d1d3_x_offset("hd1d3xo", "D1D3 x offset", 20, -10, 10);
	// T0D1,D3 y offset
	TH1F hist_d1d3_y_offset("hd1d3yo", "D1D3 y offset", 20, -10, 10);
	// T0D2,D3 x offset
	TH1F hist_d2d3_x_offset("hd2d3xo", "D2D3 x offset", 20, -10, 10);
	// T0D2,D3 y offset
	TH1F hist_d2d3_y_offset("hd2d3yo", "D2D3 y offset", 20, -10, 10);
	// detecting efficiency
	TGraph g_efficiency[3];
	// detecting efficiency multi graph
	TMultiGraph mg_efficiency;
	// energy resolution
	TH1F hist_energy_res("heres", "energy resolution", 100, -3, 3);
	// output tree
	TTree detect_tree("tree", "detected data");
	// output detect data
	DetectEvent detect;
	// setup output branches
	detect.SetupOutput(&detect_tree);

	// initialize calculators
	elc::RangeEnergyCalculator be10_calculator("10Be", "Si");
	elc::RangeEnergyCalculator he4_calculator("4He", "Si");
	elc::RangeEnergyCalculator h2_calculator("2H", "Si");
	elc::RangeEnergyCalculator* frag_calculators[2] = {
		&be10_calculator, &he4_calculator
	};

	// initialize random number generator
	TRandom3 generator(0);

	// efficiency point x
	double efficiency_x[3][100];
	// efficiency point y
	double efficiency_y[3][100];
	// counts of generated data in a specific 10Be and 14C excited energy range
	int generate_counts[3][100];
	// counts of detected data in a specific 10Be and 14C excited energy range
	int detect_counts[3][100];
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 100; ++j) {
			efficiency_x[i][j] = state_start[i] + j * 0.5;
			generate_counts[i][j] = 0;
			detect_counts[i][j] = 0;
		}
	}

	// show start
	printf("Simulating   0%%");
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

		// calculate taf energy
		double recoil_range = h2_calculator.Range(event.recoil_kinematic);
		if (recoil_range < 150.0 / cos(event.recoil_theta)) {
			detect.taf_layer = 0;
			detect.taf_lost_energy[0] = event.recoil_kinematic;
			detect.taf_lost_energy[1] = 0.0;
		} else {
			detect.taf_layer = 1;
			detect.taf_lost_energy[1] = h2_calculator.Energy(
				recoil_range - 150.0 / cos(event.recoil_theta)
			);
			detect.taf_lost_energy[0] =
				event.recoil_kinematic - detect.taf_lost_energy[1];
		}
		// detected energy
		detect.taf_energy[0] =
			detect.taf_lost_energy[0] + generator.Gaus(0.0, 0.1);
		detect.taf_energy[1] =
			detect.taf_lost_energy[1] + generator.Gaus(0.0, 1.0);

		// calculate t0 energy
		double range[2];
		double residual_energy[2];
		for (size_t i = 0; i < 2; ++i) {
			for (size_t j = 0; j < 7; ++j) {
				detect.t0_lost_energy[i][j] = 0.0;
			}
			range[i] = frag_calculators[i]->Range(event.fragment_kinematic[i]);
			residual_energy[i] = event.fragment_kinematic[i];

			for (size_t j = 0; j < 6; ++j) {
				double thick = t0_thickness[j] / cos(event.fragment_theta[i]);
				if (range[i] < thick) {
					detect.t0_layer[i] = j;
					detect.t0_lost_energy[i][j] = residual_energy[i];
					break;
				} else {
					detect.t0_lost_energy[i][j] =
						residual_energy[i]
						- frag_calculators[i]->Energy(range[i] - thick);
					range[i] -= thick;
					residual_energy[i] -= detect.t0_lost_energy[i][j];
					if (j == 5) {
						detect.t0_layer[i] = 6;
						detect.t0_lost_energy[i][6] = residual_energy[i];
					}
				}
			}
		}

		// simulate detected energy
		for (int i = 0; i < 2; ++i) {
			detect.t0_energy[i][0] =
				detect.t0_lost_energy[i][0] + generator.Gaus(0.0, 0.2);
			detect.t0_energy[i][1] =
				detect.t0_lost_energy[i][1] + generator.Gaus(0.0, 0.8);
			detect.t0_energy[i][2] =
				detect.t0_lost_energy[i][2] + generator.Gaus(0.0, 0.2);
			detect.t0_energy[i][3] =
				detect.t0_lost_energy[i][3] + generator.Gaus(0.0, 0.1);
			detect.t0_energy[i][4] =
				detect.t0_lost_energy[i][4] + generator.Gaus(0.0, 0.1);
			detect.t0_energy[i][5] =
				detect.t0_lost_energy[i][5] + generator.Gaus(0.0, 0.1);
			detect.t0_energy[i][6] =
				detect.t0_lost_energy[i][6] + generator.Gaus(0.0, 1.0);
		}

		// calculate position
		// TAF position
		// tafd strip width
		double taf_ring_width = (170.5 - 68.0) / 16.0;
		// taf ring strip
		int taf_ring_strip = int((event.rr - 68.0) / taf_ring_width);
		// taf r
		detect.tafr = (taf_ring_strip + 0.5) * taf_ring_width + 68.0;
		// a virtual taf index
		int taf_index = int(((event.recoil_phi / pi) + 1) * 3.0);
		// effective taf phi, normalized to 0-60 degree
		double effective_taf_phi =
			(event.recoil_phi / pi + 1.0) * 180.0 - taf_index * 60.0;
		// cover range of one phi strip
		double taf_phi_width = 55.2 / 8.0;
		// taf phi strip
		int taf_phi_strip = int((effective_taf_phi - 2.4) / taf_phi_width);
		// get the phi range
		double taf_phi = (
			(taf_phi_strip + 0.5) * taf_phi_width
			+ 2.4 + taf_index * 60.0 - 180.0
		) / 180.0 * pi;
		detect.tafz = event.rz;
		detect.tafx = detect.tafr * cos(taf_phi);
		detect.tafy = detect.tafr * sin(taf_phi);
		if (!taf_strips_affect) {
			detect.tafx = event.rx;
			detect.tafy = event.ry;
		}

		// T0 position
		for (size_t i = 0; i < 2; ++i) {
			// T0D1 position
			int t0d1_x_strip = int(event.fragment_x[i] + 32.0);
			detect.t0x[i][0] = double(t0d1_x_strip) - 31.5;
			int t0d1_y_strip = int(event.fragment_y[i] + 32.0);
			detect.t0y[i][0] = double(t0d1_y_strip) - 31.5;
			if (!t0_strips_affect) {
				detect.t0x[i][0] = event.fragment_x[i];
				detect.t0y[i][0] = event.fragment_y[i];
			}
			detect.t0z[i][0] = event.fragment_z[i];
			detect.t0r[i][0] = sqrt(
				pow(detect.t0x[i][0], 2.0) + pow(detect.t0y[i][0], 2.0)
			);

			// T0D2, D3 position
			for (size_t j = 1; j < 3; ++j) {
				double z = j == 1 ? 111.76 : 123.52;
				double r = z * tan(event.fragment_theta[i]);
				double x = r * cos(event.fragment_phi[i]) + event.target_x;
				double y = r * sin(event.fragment_phi[i]) + event.target_y;
				int xstrip = int((x+32.0)/2.0);
				int ystrip = int((y+32.0)/2.0);
				detect.t0x[i][j] = double(xstrip) * 2.0 - 31.0;
				detect.t0y[i][j] = double(ystrip) * 2.0 - 31.0;
				detect.t0z[i][j] = z;
				detect.t0r[i][j] = sqrt(
					pow(detect.t0x[i][j], 2.0) + pow(detect.t0y[i][j], 2.0)
				);
			}
		}

		// PPAC position
		for (size_t i = 0; i < 3; ++i) {
			detect.ppacx[i] =
				ppac_xz[i] * tan(event.beam_theta) * cos(event.beam_phi)
				+ event.target_x;
			// simulated strips
			detect.ppacx[i] = int(detect.ppacx[i] + 50.5) - 50.0;
			detect.ppacy[i] =
				ppac_yz[i] * tan(event.beam_theta) * sin(event.beam_phi)
				+ event.target_y;
			// simulated strips
			detect.ppacy[i] = int(detect.ppacy[i] + 50.5) - 50.0;
		}
		// fit and get reaction point
		double kx, ky;
		SimpleFit(ppac_xz, detect.ppacx, kx, detect.tx);
		SimpleFit(ppac_yz, detect.ppacy, ky, detect.ty);
		if (!ppac_strips_affect) {
			detect.tx = event.target_x;
			detect.ty = event.target_y;
		}

		// check valid
		detect.valid = 0;
		if (
			event.rr > 68.0 && event.rr < 170.5
			&& effective_taf_phi > 2.4 && effective_taf_phi < 57.6
			&& detect.taf_layer == 1
			&& detect.taf_lost_energy[1] > 6
		) {
			detect.valid |= 0x1;
		}
		// check T0 geometry range
		for (size_t i = 0; i < 2; ++i) {
			double l = 100.0 + detect.t0_layer[i] * 11.76;
			double x =
				l * tan(event.fragment_theta[i]) * cos(event.fragment_phi[i]);
			double y =
				l * tan(event.fragment_theta[i]) * cos(event.fragment_phi[i]);
			if (
				detect.t0_layer[i] > 0
				&& x > -32.0 && x < 32.0
				&& y > -32.0 && y < 32.0
			) {
				detect.valid |= 0x1 << (i+1);
			}
		}
		// check T0 hit same strip
		// check T0D1, T0D2
		if (
			detect.t0_layer[0] > 0 && detect.t0_layer[1] > 0
			&& (
				detect.t0x[0][0] == detect.t0x[1][0]
				|| detect.t0y[0][0] == detect.t0y[1][0]
				|| detect.t0x[0][1] == detect.t0x[1][1]
				|| detect.t0y[0][1] == detect.t0y[1][1]
			)
		) {
			detect.valid &= 0x1;
		}
		// check T0D3
		if (
			detect.t0_layer[0] > 1 && detect.t0_layer[1] > 1
			&& (
				detect.t0x[0][2] == detect.t0x[1][2]
				|| detect.t0y[0][2] == detect.t0y[1][2]
			)
		) {
			detect.valid &= 0x1;
		}

		// rebuild kinematic energy
		double taf_kinematic;
		if (taf_energy_affect) {
			taf_kinematic = detect.taf_energy[0];
			taf_kinematic += detect.taf_layer == 1 ? detect.taf_energy[1] : 0.0;
		} else {
			taf_kinematic = event.recoil_kinematic;
		}
		double t0_kinematic[2] = {0.0, 0.0};
		if (t0_energy_affect) {
			for (size_t i = 0; i < 2; ++i) {
				for (int j = 0; j < detect.t0_layer[i]+1; ++j) {
					t0_kinematic[i] += detect.t0_energy[i][j];
				}
			}
		} else {
			t0_kinematic[0] = event.fragment_kinematic[0];
			t0_kinematic[1] = event.fragment_kinematic[1];
		}

		// rebuild Q value spectrum
		ROOT::Math::XYZVector fp[2];
		for (size_t i = 0; i < 2; ++i) {
			fp[i] = ROOT::Math::XYZVector(
				detect.t0x[i][0] - detect.tx,
				detect.t0y[i][0] - detect.ty,
				detect.t0z[i][0]
			);
			double mass = i == 0 ? be10_mass : he4_mass;
			double momentum = sqrt(
				pow(t0_kinematic[i], 2.0) + 2.0 * t0_kinematic[i] * mass
			);
			fp[i] = fp[i].Unit() * momentum;
		}
		double recoil_momentum = sqrt(
			pow(taf_kinematic, 2.0) + 2.0 * taf_kinematic * h2_mass
		);
		ROOT::Math::XYZVector rp(
			detect.tafx - detect.tx,
			detect.tafy - detect.ty,
			detect.tafz
		);
		rp = rp.Unit() * recoil_momentum;
		ROOT::Math::XYZVector bp = fp[0] + fp[1] + rp;
		double calculated_beam_kinematic =
			sqrt(bp.Dot(bp) + pow(c14_mass, 2.0)) - c14_mass;
		// Q value
		detect.q = t0_kinematic[0] + t0_kinematic[1]
			+ taf_kinematic - calculated_beam_kinematic;

		// rebuild 14C excited energy
		// rebuild state
		int rebuild_state = -1;
		double rebuild_be_excited = 0.0;
		if (detect.q > -13 && detect.q < -10) {
			rebuild_state = 0;
			rebuild_be_excited = 0.0;
		} else if (detect.q > -16 && detect.q < -14) {
			rebuild_state = 1;
			rebuild_be_excited = 3.368;
		} else if (detect.q > -19 && detect.q < -16.7) {
			rebuild_state = 2;
			rebuild_be_excited = 6.179;
		}
		// calcualted 14C momentum vecotr
		ROOT::Math::XYZVector cbp = fp[0] + fp[1];
		// excited 14C momentum
		double c_momentum = cbp.R();
		// excited 14C total energy
		double c_energy =
			(t0_kinematic[0] + be10_mass + rebuild_be_excited)
			+ (t0_kinematic[1] + he4_mass);
		// excited 14C mass
		double excited_c_mass = sqrt(
			pow(c_energy, 2.0) - pow(c_momentum, 2.0)
		);
		// excited energy of 14C
		double excited_14c = excited_c_mass - c14_mass;
		if (detect.valid == 7 && rebuild_state != -1) {
			hist_energy_res.Fill(excited_14c - event.beam_excited_energy);
		}
		// fill to tree
		detect_tree.Fill();

		// increase counts
		int state = 0;
		if (event.fragment_excited_energy == 3.368) {
			state = 1;
		} else if (event.fragment_excited_energy == 6.179) {
			state = 2;
		}
		int index = int(round(
			(event.beam_excited_energy - state_start[state]) / 0.5
		));
		if (index >= 0 && index < 100) {
			generate_counts[state][index]++;
			if (detect.valid == 7 && rebuild_state != -1) detect_counts[state][index]++;
		}

		// fill DSSD offset histograms
		for (int i = 0; i < 2; ++i) {
			hist_d1d2_x_offset.Fill(detect.t0x[i][0] - detect.t0x[i][1]);
			hist_d1d2_y_offset.Fill(detect.t0y[i][0] - detect.t0y[i][1]);
			hist_d1d3_x_offset.Fill(detect.t0x[i][0] - detect.t0x[i][2]);
			hist_d1d3_y_offset.Fill(detect.t0y[i][0] - detect.t0y[i][2]);
			hist_d2d3_x_offset.Fill(detect.t0x[i][1] - detect.t0x[i][2]);
			hist_d2d3_y_offset.Fill(detect.t0y[i][1] - detect.t0y[i][2]);
		}

		// rebuild fake Q value
		int rc = 12;
		int pc = 12;
		for (int i = 0; i <= rc; ++i) {
			double fake_tafr =
				detect.tafr
				+ (i - rc/2) / double(rc) * taf_ring_width;
			for (int j = 0; j <= pc; ++j) {
				double fake_taf_phi = taf_phi
					+ (j - pc/2) / double(pc) * taf_phi_width / 180.0 * pi;
				double fake_tafx = fake_tafr * cos(fake_taf_phi);
				double fake_tafy = fake_tafr * sin(fake_taf_phi);
				// fake_tafx = tafr * cos(taf_phi);
				// fake_tafy = tafr * sin(taf_phi);
				ROOT::Math::XYZVector fake_rp(
					fake_tafx-detect.tx, fake_tafy-detect.ty, detect.tafz
				);
				fake_rp = fake_rp.Unit() * recoil_momentum;
				ROOT::Math::XYZVector fake_bp = fp[0] + fp[1] + fake_rp;
				double fake_bk =
					sqrt(fake_bp.Dot(fake_bp) + pow(c14_mass, 2.0)) - c14_mass;
				double fake_q = t0_kinematic[0] + t0_kinematic[1]
					+ taf_kinematic - fake_bk;
				if (detect.valid==7) hist_fake_q.Fill(fake_q);
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// calculate efficiency
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 100; ++j) {
			if (generate_counts[i][j] == 0) continue;
			efficiency_y[i][j] =
				double(detect_counts[i][j]) / generate_counts[i][j];
			g_efficiency[i].AddPoint(
				efficiency_x[i][j], efficiency_y[i][j]
			);
		}
	}
	mg_efficiency.Add(g_efficiency);
	mg_efficiency.Add(g_efficiency+1);
	mg_efficiency.Add(g_efficiency+2);

	// fit energy resolution
	TF1 fit_res("f1", "gaus", -2, 2);
	hist_energy_res.Fit(&fit_res, "R+");
	std::cout << "Energy resolution " <<
		fit_res.GetParameter(2) * 2.0 * sqrt(2.0 * log(2.0)) << "\n";


	detect_file.cd();
	// save histogram
	hist_fake_q.Write();
	hist_d1d2_x_offset.Write();
	hist_d1d2_y_offset.Write();
	hist_d1d3_x_offset.Write();
	hist_d1d3_y_offset.Write();
	hist_d2d3_x_offset.Write();
	hist_d2d3_y_offset.Write();
	for (int i = 0; i < 3; ++i) {
		g_efficiency[i].Write(TString::Format(
			"gef%d", i
		));
	}
	mg_efficiency.Write("mgef");
	hist_energy_res.Write();
	// save tree
	detect_tree.Write();
	// close file
	detect_file.Close();
	generate_file.Close();

	return 0;
}