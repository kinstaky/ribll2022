#include <iostream>
#include <vector>

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <Math/Vector3D.h>

#include "include/event/channel_v2_event.h"
#include "include/event/generate_event.h"
#include "include/calculator/target_energy_calculator.h"
#include "include/event/particle_event.h"
#include "include/ppac_track.h"

constexpr double TARGET_THICK = 9.53;
constexpr double POS_X_ERROR = 0.0;
constexpr double POS_Y_ERROR = 0.0;
constexpr double POS_Z_ERROR = 2.0;
constexpr double ENERGY_ERROR = 0.0;


constexpr double exbins = 200;
constexpr double exborder = 2;

using namespace ribll;

int main() {
	// input file name
	TString channel_file_name = TString::Format(
		"%s%sC14-10Be-4He-2H-v2-sim.root",
		kGenerateDataPath,
		kChannelDir
	);
	// channel file
	TFile channel_file(channel_file_name, "read");
	// input tree
	TTree *channel_tree = (TTree*)channel_file.Get("tree");
	if (!channel_tree) {
		std::cerr << "Error: Get tree from "
			<< channel_file_name << " failed.\n";
        return -1;
	}
	// setup input branches
	ChannelV2Event channel;
    channel.SetupInput(channel_tree);

	// generate file name
	TString generate_file_name = TString::Format(
        "%s%sgenerate-0004.root",
        kGenerateDataPath,
        kSimulateDir
    );
    // generate file
    TFile generate_file(generate_file_name, "read");
    // input tree
    TTree *generate_tree = (TTree*)generate_file.Get("tree");
	if (!generate_tree) {
		std::cerr << "Error: Get tree from "
			<< generate_file_name << " failed.\n";
		return -1;
	}
	GenerateEvent generate;
	// setup input branches
	generate.SetupInput(generate_tree);


    // output file name
    TString output_file_name = TString::Format(
        "%s%ssystematic-error.root",
        kGenerateDataPath,
        kSimulateDir
    );
    // output file
    TFile opf(output_file_name, "recreate");
	// excitation energy difference
	TH1F hdex[3] = {
		TH1F(
			"hdex0", "excitation energy difference of g.s.",
			exbins, -exborder, exborder
		),
		TH1F(
			"hdex1", "excitation energy difference of 2_1^+",
			exbins, -exborder, exborder
		),
		TH1F(
			"hdex2", "excitation energy difference of 6MeV",
			exbins, -exborder, exborder
		)
	};
	TH1F hdmex[3] = {
		TH1F(
			"hdmex0", "18MeV excitation energy difference of g.s.",
			exbins, -exborder, exborder
		),
		TH1F(
			"hdmex1", "21.5MeV excitation energy difference of 2_1^+",
			exbins, -exborder, exborder
		),
		TH1F(
			"hdmex2", "24MeV excitation energy difference of 6MeV",
			exbins, -exborder, exborder
		)
	};
	TH2F hdex2[3] = {
		TH2F(
			"hdexe0", "excitation energy difference of g.s. VS ex energy",
			100, 12, 32,
			exbins, -exborder, exborder
		),
		TH2F(
			"hdexe1", "excitation energy difference of 2_1^+ VS ex energy",
			100, 15, 35,
			exbins, -exborder, exborder
		),
		TH2F(
			"hdexe2", "excitation energy difference of 6MeV VS ex energy",
			100, 18, 38,
			exbins, -exborder, exborder
		)
	};

	elc::TargetEnergyCalculator be10_target("10Be", "CD2", TARGET_THICK);
	elc::TargetEnergyCalculator he4_target("4He", "CD2", TARGET_THICK);


	long long entries = channel_tree->GetEntries();
	long long entry100 = entries / 100 + 1;
	printf("Processing   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < channel_tree->GetEntriesFast(); ++entry) {
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		channel_tree->GetEntry(entry);
		if (channel.valid != 0) continue;
		generate_tree->GetEntry(channel.entry);


		double stx = channel.tx;
		double sty = channel.ty;
		double sbex = channel.fragment_x[0];
		double sbey = channel.fragment_y[0];
		double sbez = 100.0;
		double shex = channel.fragment_x[1];
		double shey = channel.fragment_y[1];
		double shez = 100.0;
		double sbek = channel.fragment_kinetic[0];
		double shek = channel.fragment_kinetic[1];
		double mass_excited_10be = mass_10be + generate.fragment_excited_energy;

		// consider systematic error of T0 position
		sbex += POS_X_ERROR;
		shex += POS_X_ERROR;
		sbey += POS_Y_ERROR;
		shey += POS_Y_ERROR;
		sbez += POS_Z_ERROR;
		shez += POS_Z_ERROR;

		// get directions
		ROOT::Math::XYZVector be_direction(
			sbex - stx, sbey - sty, sbez
		);
		be_direction = be_direction.Unit();

		ROOT::Math::XYZVector he_direction(
			shex - stx, shey - sty, shez
		);
		he_direction = he_direction.Unit();


		// consider systematic error of T0 energy
		sbek += ENERGY_ERROR;
		shek += ENERGY_ERROR;

		// consider target thick
		sbek = be10_target.Energy(
			-0.5 / cos(be_direction.Theta()),
			sbek
		);
		shek = he4_target.Energy(
            -0.5 / cos(he_direction.Theta()),
            shek
        );

		ROOT::Math::XYZVector bep = be_direction * MomentumFromKinetic(
			mass_excited_10be, sbek
		);
		ROOT::Math::XYZVector hep = he_direction * MomentumFromKinetic(
			mass_4he, shek
		);

		// calculate excited energy
		ROOT::Math::XYZVector xcp = bep + hep;
		double xc_energy = sbek + mass_excited_10be + shek + mass_4he;
		double xc_mass = sqrt(pow(xc_energy, 2.0) - xcp.Mag2());
		double xce = xc_mass - mass_14c;

		hdex[generate.fragment_state].Fill(xce-generate.beam_excited_energy);
		if (generate.fragment_state == 0) {
			if (
				generate.beam_excited_energy > 17.5
				&& generate.beam_excited_energy < 18.5
			) {
				hdmex[0].Fill(xce-generate.beam_excited_energy);
			}
		} else if (generate.fragment_state == 1) {
			if (
				generate.beam_excited_energy > 21.0
				&& generate.beam_excited_energy < 22.0
			) {
				hdmex[1].Fill(xce-generate.beam_excited_energy);
			}
		} else if (generate.fragment_state == 2) {
			if (
				generate.beam_excited_energy > 23.5
				&& generate.beam_excited_energy < 24.5
			) {
				hdmex[2].Fill(xce-generate.beam_excited_energy);
			}
		}
		hdex2[generate.fragment_state].Fill(
			generate.beam_excited_energy, xce-generate.beam_excited_energy
		);
	}
	printf("\b\b\b\b100%%\n");

	// fit
	for (int i = 0; i < 3; ++i) {
		TF1 *f1 = new TF1(TString::Format("fx%d", i), "gaus", -exborder, exborder);
		f1->SetNpx(10000);
		hdex[i].Fit(f1, "QR+");
		std::cout << "dEx " << i << ": "
			<< f1->GetParameter(1)*1e3 << ", " << f1->GetParameter(2)*1e3 << "\n";
	}

	for (int i = 0; i < 3; ++i) {
		TF1 *f1 = new TF1(TString::Format("fmx%d", i), "gaus", -exborder, exborder);
		f1->SetNpx(10000);
		hdmex[i].Fit(f1, "QR+");
		std::cout << "dmEx " << i << ": "
			<< f1->GetParameter(1)*1e3 << ", " << f1->GetParameter(2)*1e3 << "\n";
	}

	opf.cd();
	for (size_t i = 0; i < 3; ++i) hdex[i].Write();
	for (size_t i = 0; i < 3; ++i) hdmex[i].Write();
	for (size_t i = 0; i < 3; ++i) hdex2[i].Write();
	opf.Close();

	generate_file.Close();
	channel_file.Close();

	return 0;
}