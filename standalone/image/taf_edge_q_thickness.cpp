#include <iostream>

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <Math/Vector3D.h>

#include "include/event/channel_v2_event.h"
#include "include/calculator/csi_energy_calculator.h"

using namespace ribll;

int main() {
	// input file name
	TString input_file_name = TString::Format(
		"%s%sC14-10Be-4He-2H-v2.root",
		kGenerateDataPath,
		kChannelDir
	);
	// input file
	TFile ipf(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// input event
	ChannelV2Event channel;
	// setup input branches
	channel.SetupInput(ipt);

	std::vector<double> thickness;
	for (double t = 140.0; t <= 180.0; t += 2.0) thickness.push_back(t);

	// histograms
	std::vector<TH1F> hist_q_all_strips[6];
	std::vector<TH1F> hist_q_single_strip[6][3];
	for (const auto &t : thickness) {
		for (int i = 0; i < 6; ++i) {
			hist_q_all_strips[i].emplace_back(
				TString::Format("hq%dt%d", i, int(t)),
				TString::Format(
					"Q spectrum of TAF %d at thickness %d um",
					i, int(t)
				),
				30, -23, -8
			);
			for (int j = 0; j < 3; ++j) {
				hist_q_single_strip[i][j].emplace_back(
					TString::Format("hq%ds%dt%d", i, j+13, int(t)),
					TString::Format(
						"Q spectrum of TAF %d strip %d at thickness %d",
						i, j+13, int(t)
					),
					30, -23, -8
				);
			}
		}
	}


	// get mass
	constexpr double beam_mass = mass_14c;
	constexpr double recoil_mass = mass_2h;

	// calculators
	elc::CsiEnergyCalculator recoil_calculator("2H");

	// loop
	for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
		ipt->GetEntry(entry);

		if (channel.valid != 0) continue;
		if (channel.tafd_front_strip < 13) continue;

		int &index = channel.taf_index;
		int &fs = channel.tafd_front_strip;

		// get recoil direction
		ROOT::Math::XYZVector recoil_direction(
			channel.recoil_x - channel.tx,
			channel.recoil_y - channel.ty,
			channel.recoil_z
		);
		recoil_direction = recoil_direction.Unit();

		// get fragment momentum vector
		ROOT::Math::XYZVector fragment_p[2];
		for (int i = 0; i < 2; ++i) {
			fragment_p[i] = ROOT::Math::XYZVector(
				channel.fragment_x[i] - channel.tx,
				channel.fragment_y[i] - channel.ty,
				channel.fragment_z[i]
			);
			fragment_p[i] =
				fragment_p[i].Unit() * channel.fragment_momentum[i];
		}

		for (size_t i = 0; i < thickness.size(); ++i) {
			// get recoil energy
			double recoil_kinetic = recoil_calculator.Energy(
				recoil_direction.Theta(),
				channel.tafd_energy,
				thickness[i]
			);
			double recoil_momentum = MomentumFromKinetic(
				recoil_mass, recoil_kinetic
			);
			ROOT::Math::XYZVector recoil_p = recoil_direction * recoil_momentum;

			// get beam momentum and energy
			ROOT::Math::XYZVector beam_p =
				fragment_p[0] + fragment_p[1] + recoil_p;
			double beam_kinetic = KineticFromMomentum(
				beam_mass, beam_p.R()
			);

			double q = channel.fragment_kinetic[0] + channel.fragment_kinetic[1]
				+ recoil_kinetic - beam_kinetic;

			hist_q_all_strips[index][i].Fill(q);
			hist_q_single_strip[index][fs-13][i].Fill(q);
		}
	}

	TCanvas *c1 = new TCanvas("c1", "c1", 1920, 1080);
	c1->cd();

	TString all_strips_pdf_file_name = TString::Format(
		"%s%staf-edge-q-thickness-all.pdf",
		kGenerateDataPath,
		kImageDir
	);
	c1->Print(all_strips_pdf_file_name+"[");
	for (size_t i = 0; i < 6; ++i) {
		for (auto &hist : hist_q_all_strips[i]) {
			hist.Draw();
			c1->Print(all_strips_pdf_file_name);
		}
	}
	c1->Print(all_strips_pdf_file_name+"]");

	TString single_strip_pdf_file_name = TString::Format(
		"%s%staf-edge-q-thickness-single.pdf",
		kGenerateDataPath,
		kImageDir
	);
	c1->Print(single_strip_pdf_file_name+"[");
	for (size_t i = 0; i < 6; ++i) {
		for (size_t j = 0; j < 3; ++j) {
			for (auto &hist : hist_q_single_strip[i][j]) {
				hist.Draw();
				c1->Print(single_strip_pdf_file_name);
			}
		}
	}
	c1->Print(single_strip_pdf_file_name+"]");


	return 0;
}