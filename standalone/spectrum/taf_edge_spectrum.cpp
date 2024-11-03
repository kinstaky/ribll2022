#include <iostream>

#include <TH1F.h>
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TChain.h>
#include <Math/Vector3D.h>

#include "include/event/ta_event.h"
#include "include/event/particle_event.h"
#include "include/ppac_track.h"
#include "include/calculator/csi_energy_calculator.h"

using namespace ribll;

int main() {
	// T0 chain
	TChain t0_chain("t0", "t0");
	// TAF chain
	TChain taf_chain[6] = {
		TChain{"taf0", "taf0"},
		TChain{"taf1", "taf1"},
		TChain{"taf2", "taf2"},
		TChain{"taf3", "taf3"},
		TChain{"taf4", "taf4"},
		TChain{"taf5", "taf5"}
	};
	// PPAC chain
	TChain ppac_chain("ppac", "ppac");

	for (int run = 618; run <= 716; ++run) {
		if (run == 628) continue;
		if (run > 652 && run < 675) continue;
		if (run > 716 && run < 739) continue;
		// T0 chain
		t0_chain.AddFile(TString::Format(
			"%s%st0-particle-ta-%04d.root/tree",
			kGenerateDataPath,
			kParticleDir,
			run
		));
		// TAF chain
		for (int i = 0; i < 6; ++i) {
			taf_chain[i].AddFile(TString::Format(
				"%s%staf%d-telescope-ta-%04d.root/tree",
				kGenerateDataPath,
				kTelescopeDir,
				i,
				run
			));
		}
		// PPAC chain
		ppac_chain.AddFile(TString::Format(
			"%s%sxppac-particle-ta-%04d.root/tree",
			kGenerateDataPath,
			kParticleDir,
			run
		));
	}
	// add friends
	for (int i = 0; i < 6; ++i) {
		t0_chain.AddFriend(taf_chain+i, TString::Format("taf%d", i));
	}
	t0_chain.AddFriend(&ppac_chain, "ppac");
	// input T0 events
	ParticleEvent t0_event;
	// input TAF events
	TaEvent taf_event[6];
	// input PPAC event
	ParticleEvent ppac_event;
	unsigned short ppac_xflag, ppac_yflag;
	// setup input branches
	t0_event.SetupInput(&t0_chain);
	for (int i = 0; i < 6; ++i) {
		taf_event[i].SetupInput(
			&t0_chain, TString::Format("taf%d.", i).Data()
		);
	}
	ppac_event.SetupInput(&t0_chain, "ppac.");
	t0_chain.SetBranchAddress("ppac.xflag", &ppac_xflag);
	t0_chain.SetBranchAddress("ppac.yflag", &ppac_yflag);

	// output file name
	TString output_file_name = TString::Format(
		"%s%staf_edge_spectrum.root",
		kGenerateDataPath,
		kSpectrumDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "TAF edge strip spectrum");
	// output data
	int front_strip, taf_index;
	double tafde, csie, q;
	double bek, hek, dk, ck;
	bool csi_valid;
	opt.Branch("taf_index", &taf_index, "tafi/I");
	opt.Branch("front_strip", &front_strip, "fs/I");
	opt.Branch("csi_valid", &csi_valid, "cv/O");
	opt.Branch("tafde", &tafde, "tafde/D");
	opt.Branch("csie", &csie, "csie/D");
	opt.Branch("be_kinetic", &bek, "bek/D");
	opt.Branch("he_kinetic", &hek, "hek/D");
	opt.Branch("d_kinetic", &dk, "dk/D");
	opt.Branch("c_kinetic", &ck, "ck/D");
	opt.Branch("q", &q, "q/D");

	elc::CsiEnergyCalculator h2_calculator("2H");

	int conflict_deutron = 0;

	// total number of entries
	long long entries = t0_chain.GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Filling TAF pid   0%%");
	fflush(stdout);
	// loop to fill pid histogram
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		// get event
		t0_chain.GetEntry(entry);

		// check Be and He
		int he_index = -1;
		int be_index = -1;
		for (int i = 0; i < t0_event.num; ++i) {
			if (t0_event.charge[i] == 2 && t0_event.mass[i] == 4) he_index = i;
			if (t0_event.charge[i] == 4 && t0_event.mass[i] == 10) be_index = i;
		}
		if (he_index < 0 || be_index < 0) continue;
		// check possible deutron
		int deutron_number = 0;
		taf_index = -1;
		for (int i = 0; i < 6; ++i) {
			if (taf_event[i].num == 0) continue;
			taf_index = i;
			++deutron_number;
			if (taf_event[i].flag[0] == 1) csi_valid = false;
			else csi_valid = true;
		}
		if (taf_index < 0) continue;
		if (deutron_number > 1) {
			++conflict_deutron;
			continue;
		}
		front_strip = taf_event[taf_index].front_strip[0];
		// check PPAC
		if (ppac_event.num != 4) continue;

		// PPAC
		// get using PPAC xz and yz
		double using_ppac_xz[3] = {ppac_xz[0], ppac_xz[1], ppac_xz[2]};
		double using_ppac_yz[3] = {ppac_yz[0], ppac_yz[1], ppac_yz[2]};
		// if (run >= ppac_change_run) {
		// 	using_ppac_xz[0] = all_ppac_xz[1];
		// 	using_ppac_yz[0] = all_ppac_yz[1];
		// }
		// get using ppac correct
		double using_ppac_correct[2][3] = {
			{ppac_correct[0][0], ppac_correct[0][1], ppac_correct[0][2]},
			{ppac_correct[1][0], ppac_correct[1][1], ppac_correct[1][2]}
		};
		// if (run >= ppac_change_run) {
		// 	using_ppac_correct[0][0] = all_ppac_correct[0][1];
		// 	using_ppac_correct[1][0] = all_ppac_correct[1][1];
		// }
		double xppac_x[3], xppac_y[3];
		for (int i = 0; i < 3; ++i) {
			xppac_x[i] = ppac_event.x[i] - using_ppac_correct[0][i];
			xppac_y[i] = ppac_event.y[i] - using_ppac_correct[1][i];
		}

		double xk, yk, tx, ty;
		TrackMultiplePpac(ppac_xflag, using_ppac_xz, xppac_x, xk, tx);
		TrackMultiplePpac(ppac_yflag, using_ppac_yz, xppac_y, yk, ty);

		// 10Be momentum
		double be_momentum = MomentumFromKinetic(mass_10be, t0_event.energy[be_index]);
		// 10Be momentum vector
		ROOT::Math::XYZVector p_be(
			t0_event.x[be_index] - tx,
			t0_event.y[be_index] - ty,
			100.0
		);
		p_be = p_be.Unit() * be_momentum;

		// 4He momentum
		double he_momentum = MomentumFromKinetic(mass_4he, t0_event.energy[he_index]);
		// 4He momentum vector
		ROOT::Math::XYZVector p_he(
			t0_event.x[he_index] - tx,
			t0_event.y[he_index] - ty,
			100.0
		);
		p_he = p_he.Unit() * he_momentum;

		// 2H momentum vector
		ROOT::Math::Polar3DVector d_position(
			taf_event[taf_index].radius[0],
			taf_event[taf_index].theta[0],
			taf_event[taf_index].phi[0]
		);
		ROOT::Math::XYZVector p_d(
			d_position.X() - tx,
			d_position.Y() - ty,
			135.0
		);
		// TAFD energy
		tafde = taf_event[taf_index].energy[0][0];
		// TAF CsI energy
		csie = h2_calculator.Energy(
			p_d.Theta(), tafde, tafd_thickness[taf_index]
		);
		double taf_energy = tafde + csie;
		// 2H momentum
		double d_momentum = MomentumFromKinetic(mass_2h, taf_energy);
		p_d = p_d.Unit() * d_momentum;

		// beam 14C momentum vector
		ROOT::Math::XYZVector p_beam = p_be + p_he + p_d;

		// 14C momentum
		double c14_momentum = p_beam.R();
		// 14C kinematic energy
		double c14_kinetic =
			sqrt(pow(c14_momentum, 2.0) + pow(mass_14c, 2.0)) - mass_14c;

		// three-fold Q value
		q = t0_event.energy[be_index] + t0_event.energy[he_index]
			+ taf_energy - c14_kinetic;

		bek = t0_event.energy[be_index];
		hek = t0_event.energy[he_index];
		dk = taf_energy;
		ck = c14_kinetic;

		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	std::cout << "Multiple deutron " << conflict_deutron << " / "
		<< opt.GetEntries() << " "
		<< double(conflict_deutron) / opt.GetEntries() << "\n";


	// save
	opf.cd();
	opt.Write();
	opf.Close();

	return 0;
}