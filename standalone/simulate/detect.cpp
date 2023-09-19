#include <cmath>
#include <iostream>

#include <TFile.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/generate_event.h"
#include "include/calculator/range_energy_calculator.h"

using namespace ribll;

constexpr double pi = 3.1415926;
constexpr double u = 931.494;
constexpr double c14_mass = 13.9999505089 * u;
constexpr double be10_mass = 10.0113403769 * u;
constexpr double he4_mass = 4.0015060943 * u;
constexpr double h2_mass = 2.0135531980 * u;

constexpr double ppac_xz[3] = {-695.2, -454.2, -275.2};
constexpr double ppac_yz[3] = {-689.2, -448.2, -269.2};

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
	TH1F hist_fake_q("hfq", "q value", 1000, -20, -10);
	// output file
	TFile detect_file(detect_file_name, "recreate");
	// output tree
	TTree detect_tree("tree", "detected data");
	// output detect data
	int taf_layer;
	double taf_lost_energy[2];
	double taf_energy[2];
	int t0_layer[2];
	double t0_lost_energy[2][7];
	double t0_energy[2][7];
	double tafx;
	double tafy;
	double tafz;
	double tafr;
	double t0x[2];
	double t0y[2];
	double t0z[2];
	double t0r[2];
	double ppacx[3];
	double ppacy[3];
	double tx;
	double ty;
	int valid;
	double q_value;
	// setup input branches
	detect_tree.Branch("taf_layer", &taf_layer, "taflayer/I");
	detect_tree.Branch("taf_lost_energy", taf_lost_energy, "tafle[2]/D");
	detect_tree.Branch("taf_energy", taf_energy, "tafe[2]/D");
	detect_tree.Branch("t0_layer", t0_layer, "t0layer[2]/I");
	detect_tree.Branch("t0_lost_energy", t0_lost_energy, "t0le[2][7]/D");
	detect_tree.Branch("t0_energy", t0_energy, "t0e[2][7]/D");
	detect_tree.Branch("tafx", &tafx, "tafx/D");
	detect_tree.Branch("tafy", &tafy, "tafy/D");
	detect_tree.Branch("tafz", &tafz, "tafz/D");
	detect_tree.Branch("tafr", &tafr, "tafr/D");
	detect_tree.Branch("t0x", t0x, "t0x[2]/D");
	detect_tree.Branch("t0y", t0y, "t0y[2]/D");
	detect_tree.Branch("t0z", t0z, "t0z[2]/D");
	detect_tree.Branch("t0r", t0r, "t0r[2]/D");
	detect_tree.Branch("ppacx", ppacx, "ppacx[3]/D");
	detect_tree.Branch("ppacy", ppacy, "ppacy[3]/D");
	detect_tree.Branch("tx", &tx, "tx/D");
	detect_tree.Branch("ty", &ty, "ty/D");
	detect_tree.Branch("valid", &valid, "valid/I");
	detect_tree.Branch("q", &q_value, "q/D");

	// initialize calculators
	elc::RangeEnergyCalculator be10_calculator("10Be", "Si");
	elc::RangeEnergyCalculator he4_calculator("4He", "Si");
	elc::RangeEnergyCalculator h2_calculator("2H", "Si");
	elc::RangeEnergyCalculator* frag_calculators[2] = {
		&be10_calculator, &he4_calculator
	};

	// initialize random number generator
	TRandom3 generator(0);

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
			taf_layer = 0;
			taf_lost_energy[0] = event.recoil_kinematic;
			taf_lost_energy[1] = 0.0;
		} else {
			taf_layer = 1;
			taf_lost_energy[1] = h2_calculator.Energy(
				recoil_range - 150.0 / cos(event.recoil_theta)
			);
			taf_lost_energy[0] = event.recoil_kinematic - taf_lost_energy[1];
		}
		// detected energy
		taf_energy[0] = taf_lost_energy[0] + generator.Gaus(0.0, 0.1);
		taf_energy[1] = taf_lost_energy[1] + generator.Gaus(0.0, 1.0);

		// calculate t0 energy
		double range[2];
		double residual_energy[2];
		for (size_t i = 0; i < 2; ++i) {
			for (size_t j = 0; j < 7; ++j) {
				t0_lost_energy[i][j] = 0.0;
			}
			range[i] = frag_calculators[i]->Range(event.fragment_kinematic[i]);
			residual_energy[i] = event.fragment_kinematic[i];

			for (size_t j = 0; j < 6; ++j) {
				double thick = t0_thickness[j] / cos(event.fragment_theta[i]);
				if (range[i] < thick) {
					t0_layer[i] = j;
					t0_lost_energy[i][j] = residual_energy[i];
					break;
				} else {
					t0_lost_energy[i][j] =
						residual_energy[i]
						- frag_calculators[i]->Energy(range[i] - thick);
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
			t0_energy[i][0] = t0_lost_energy[i][0] + generator.Gaus(0.0, 0.2);
			t0_energy[i][1] = t0_lost_energy[i][1] + generator.Gaus(0.0, 0.8);
			t0_energy[i][2] = t0_lost_energy[i][2] + generator.Gaus(0.0, 0.2);
			t0_energy[i][3] = t0_lost_energy[i][3] + generator.Gaus(0.0, 0.1);
			t0_energy[i][4] = t0_lost_energy[i][4] + generator.Gaus(0.0, 0.1);
			t0_energy[i][5] = t0_lost_energy[i][5] + generator.Gaus(0.0, 0.1);
			t0_energy[i][6] = t0_lost_energy[i][6] + generator.Gaus(0.0, 1.0);
		}

		// calculate position
		// TAF position
		// tafd strip width
		double taf_ring_width = (170.5 - 68.0) / 16.0;
		// taf ring strip
		int taf_ring_strip = int((event.rr - 68.0) / taf_ring_width);
		// taf r
		tafr = (taf_ring_strip + 0.5) * taf_ring_width + 68.0;
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
		tafz = event.rz;
		tafx = tafr * cos(taf_phi);
		tafy = tafr * sin(taf_phi);
		// tafx = event.rx;
		// tafy = event.ry;

		// T0 position
		for (size_t i = 0; i < 2; ++i) {
			int t0_x_strip = int(event.fragment_x[i] + 32.0);
			t0x[i] = double(t0_x_strip) - 31.5;
			int t0_y_strip = int(event.fragment_y[i] + 32.0);
			t0y[i] = double(t0_y_strip) - 31.5;
			// t0x[i] = event.fragment_x[i];
			// t0y[i] = event.fragment_y[i];
			t0z[i] = event.fragment_z[i];
			t0r[i] = sqrt(pow(t0x[i], 2.0) + pow(t0y[i], 2.0));
		}

		// PPAC position
		for (size_t i = 0; i < 3; ++i) {
			ppacx[i] =
				ppac_xz[i] * tan(event.beam_theta) * cos(event.beam_phi)
				+ event.target_x;
			// simulated strips
			ppacx[i] = int(ppacx[i] + 50.5) - 50.0;
			ppacy[i] =
				ppac_yz[i] * tan(event.beam_theta) * sin(event.beam_phi)
				+ event.target_y;
			// simulated strips
			ppacy[i] = int(ppacy[i] + 50.5) - 50.0;
		}
		// fit and get reaction point
		double kx, ky;
		SimpleFit(ppac_xz, ppacx, kx, tx);
		SimpleFit(ppac_yz, ppacy, ky, ty);
		// tx = event.target_x;
		// ty = event.target_y;

		// check valid
		valid = 0;
		if (
			event.rr > 68.0 && event.rr < 170.5
			&& effective_taf_phi > 2.4 && effective_taf_phi < 57.6
			&& taf_layer == 1
			&& taf_lost_energy[1] > 6
		) {
			valid |= 0x1;
		}
		for (size_t i = 0; i < 2; ++i) {
			double l = 100.0 + t0_layer[i] * 11.76;
			double x =
				l * tan(event.fragment_theta[i]) * cos(event.fragment_phi[i]);
			double y =
				l * tan(event.fragment_theta[i]) * cos(event.fragment_phi[i]);
			if (
				// t0_layer[i] > 0
				x > -32.0 && x < 32.0
				&& y > -32.0 && y < 32.0
			) {
				valid |= 0x1 << (i+1);
			}
		}

		// rebuild kinematic energy
		double taf_kinematic = taf_energy[0];
		taf_kinematic += taf_layer == 1 ? taf_energy[1] : 0.0;
		// double taf_kinematic = event.recoil_kinematic;
		// double t0_kinematic[2] = {event.fragment_kinematic[0], event.fragment_kinematic[1]};
		double t0_kinematic[2] = {0.0, 0.0};
		for (size_t i = 0; i < 2; ++i) {
			for (int j = 0; j < t0_layer[i]+1; ++j) {
				t0_kinematic[i] += t0_energy[i][j];
			}
		}

		// rebuild Q value spectrum
		ROOT::Math::XYZVector fp[2];
		for (size_t i = 0; i < 2; ++i) {
			fp[i] = ROOT::Math::XYZVector(t0x[i]-tx, t0y[i]-ty, t0z[i]);
			double mass = i == 0 ? be10_mass : he4_mass;
			double momentum = sqrt(
				pow(t0_kinematic[i], 2.0) + 2.0 * t0_kinematic[i] * mass
			);
			fp[i] = fp[i].Unit() * momentum;
		}
		double recoil_momentum = sqrt(
			pow(taf_kinematic, 2.0) + 2.0 * taf_kinematic * h2_mass
		);
		ROOT::Math::XYZVector rp(tafx-tx, tafy-ty, tafz);
		rp = rp.Unit() * recoil_momentum;
		ROOT::Math::XYZVector bp = fp[0] + fp[1] + rp;
		double calculated_beam_kinematic =
			sqrt(bp.Dot(bp) + pow(c14_mass, 2.0)) - c14_mass;
		// Q value
		q_value = t0_kinematic[0] + t0_kinematic[1]
			+ taf_kinematic - calculated_beam_kinematic;

		// fill to tree
		detect_tree.Fill();


		// rebuild fake Q value
		int rc = 12;
		int pc = 12;
		for (int i = 0; i <= rc; ++i) {
			double fake_tafr = tafr + (i - rc/2) / double(rc) * taf_ring_width;
			for (int j = 0; j <= pc; ++j) {
				double fake_taf_phi = taf_phi
					+ (j - pc/2) / double(pc) * taf_phi_width / 180.0 * pi;
				double fake_tafx = fake_tafr * cos(fake_taf_phi);
				double fake_tafy = fake_tafr * sin(fake_taf_phi);
				// fake_tafx = tafr * cos(taf_phi);
				// fake_tafy = tafr * sin(taf_phi);
				ROOT::Math::XYZVector fake_rp(
					fake_tafx-tx, fake_tafy-ty, tafz
				);
				fake_rp = fake_rp.Unit() * recoil_momentum;
				ROOT::Math::XYZVector fake_bp = fp[0] + fp[1] + fake_rp;
				double fake_bk =
					sqrt(fake_bp.Dot(fake_bp) + pow(c14_mass, 2.0)) - c14_mass;
				double fake_q = t0_kinematic[0] + t0_kinematic[1]
					+ taf_kinematic - fake_bk;
				if (valid==7) hist_fake_q.Fill(fake_q);
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");


	detect_file.cd();
	// save histogram
	hist_fake_q.Write();
	// save tree
	detect_tree.Write();
	// close file
	detect_file.Close();
	generate_file.Close();

	return 0;
}