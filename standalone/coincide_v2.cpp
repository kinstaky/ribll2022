#include <iostream>

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <Math/Vector3D.h>

#include "include/event/channel_v2_event.h"
#include "include/event/particle_event.h"
#include "include/event/ta_event.h"
#include "include/ppac_track.h"
#include "include/calculator/csi_energy_calculator.h"

using namespace ribll;

// get mass
const double beam_mass = mass_14c;
double recoil_mass = mass_2h;
const double fragment_mass[2] = {mass_10be, mass_4he};

void CoincideSingleRun(
	int run,
	bool simulation,
	TTree &opt,
	ChannelV2Event &channel,
	elc::CsiEnergyCalculator &h2_calculator,
	int &taf_conflict
) {
	// T0 file name
	TString t0_file_name = TString::Format(
		"%s%st0-particle-%sta-v2-%04d.root",
		kGenerateDataPath,
		kParticleDir,
		simulation ? "sim-" : "",
		run
	);
	// T0 file
	TFile ipf(t0_file_name, "read");
	// T0 tree
	TTree *ipt = (TTree*)ipf.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< t0_file_name << " failed.\n";
		return;
	}
	// PPAC file name
	TString ppac_file_name = TString::Format(
		"%s%sxppac-particle-%sta-%04d.root",
		kGenerateDataPath,
		kParticleDir,
		simulation ? "sim-" : "",
		run
	);
	// add friend
	ipt->AddFriend("xppac=tree", ppac_file_name);
	// TAF file names
	TString taf_file_names[6];
	for (int i = 0; i < 6; ++i) {
		taf_file_names[i].Form(
			"%s%staf%d-particle-%sta-v2-%04d.root",
			kGenerateDataPath,
			kParticleDir,
			i,
			simulation ? "sim-" : "",
			run
		);
	}
	// add friend
	for (int i = 0; i < 6; ++i) {
		ipt->AddFriend(
			TString::Format("taf%d=tree", i),
			taf_file_names[i]
		);
	}
	// TAF telescope file names
	TString taf_tele_file_names[6];
	for (int i = 0; i < 6; ++i) {
		taf_tele_file_names[i].Form(
			"%s%staf%d-telescope-%sta-%04d.root",
			kGenerateDataPath,
			kTelescopeDir,
			i,
			simulation ? "sim-" : "",
			run
		);
	}
	// add friend
	for (int i = 0; i < 6; ++i) {
		ipt->AddFriend(
			TString::Format("taf%dtele=tree", i),
			taf_tele_file_names[i]
		);
	}

	// input events
	ParticleEvent t0_event;
	bool hole[4];
	ParticleEvent xppac_event;
	unsigned short xppac_xflag, xppac_yflag;
	TaEvent taf_tele_events[6];
	ParticleEvent taf_events[6];

	// setup input branches
	t0_event.SetupInput(ipt);
	ipt->SetBranchAddress("hole", hole);
	xppac_event.SetupInput(ipt, "xppac.");
	ipt->SetBranchAddress("xppac.xflag", &xppac_xflag);
	ipt->SetBranchAddress("xppac.yflag", &xppac_yflag);
	for (int i = 0; i < 6; ++i) {
		taf_events[i].SetupInput(
			ipt, TString::Format("taf%d.", i).Data()
		);
		taf_tele_events[i].SetupInput(
			ipt, TString::Format("taf%dtele.", i).Data()
		);
	}

	// PPAC information
	// get using PPAC xz and yz
	double using_ppac_xz[3] = {ppac_xz[0], ppac_xz[1], ppac_xz[2]};
	double using_ppac_yz[3] = {ppac_yz[0], ppac_yz[1], ppac_yz[2]};
	if (run >= ppac_change_run) {
		using_ppac_xz[0] = all_ppac_xz[1];
		using_ppac_yz[0] = all_ppac_yz[1];
	}
	// get using ppac correct
	double using_ppac_correct[2][3] = {
		{ppac_correct[0][0], ppac_correct[0][1], ppac_correct[0][2]},
		{ppac_correct[1][0], ppac_correct[1][1], ppac_correct[1][2]}
	};
	if (run >= ppac_change_run) {
		using_ppac_correct[0][0] = all_ppac_correct[0][1];
		using_ppac_correct[1][0] = all_ppac_correct[1][1];
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100 + 1;
	if (simulation) {
		printf("Reading   0%%");
		fflush(stdout);
	}
	// loop
	for (long long entry = 0; entry < entries; ++entry) {
		if (simulation && entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);

		// initialize
		channel.valid = 0;
		channel.ppac_valid = false;
		channel.t0_valid = false;
		channel.tafd_edge = false;
		channel.tafcsi_valid = false;
		channel.hole = 0;

		// check Be and He
		int frag_index[2] = {-1, -1};
		for (int i = 0; i < t0_event.num; ++i) {
			if (
				t0_event.charge[i] == channel.fragment_charge[0]
				&& t0_event.mass[i] == channel.fragment_mass[0]
			) {
				frag_index[0] = i;
			}
			if (
				t0_event.charge[i] == channel.fragment_charge[1]
				&& t0_event.mass[i] == channel.fragment_mass[1]
			) {
				frag_index[1] = i;
			}
		}
		if (frag_index[0] < 0 || frag_index[1] < 0) {
			channel.t0_valid = false;
			channel.valid |= 1;
			continue;
		}
		if (hole[frag_index[0]]) channel.hole |= 1;
		if (hole[frag_index[1]]) channel.hole |= 2;
		if (channel.hole == 0) {
			channel.t0_valid = true;
		} else {
			channel.t0_valid = false;
			channel.valid |= 1;
		}

		// check recoil particle
		int taf_number = 0;
		channel.taf_index = -1;
		for (int i = 0; i < 6; ++i) {
			if (taf_events[i].num == 0) continue;
			if (taf_events[i].mass[0] != channel.recoil_mass) continue;
			channel.taf_index = i;
			++taf_number;
		}
		if (channel.taf_index < 0) {
			// check edge strip possible recoil particle
			for (int i = 0; i < 6; ++i) {
				if (taf_events[i].num == 0) continue;
				if (taf_events[i].mass[0] != 0) continue;
				if (taf_tele_events[i].front_strip[0] < 12) continue;
				channel.taf_index = i;
				++taf_number;
			}
			channel.tafd_edge = true;
		} else {
			channel.tafcsi_valid = true;
		}
		// recoil particle not found
		if (channel.taf_index < 0) {
			channel.valid |= 2;
			continue;
		}
		// found multiple recoil particle
		if (taf_number > 1) {
			++taf_conflict;
			channel.valid |= 2;
			continue;
		}
		channel.csi_index =
			channel.taf_index*2 + taf_events[channel.taf_index].index[0];

		// check PPAC
		channel.ppac_xflag = xppac_xflag;
		channel.ppac_yflag = xppac_yflag;
		channel.ppac_xnum = channel.ppac_ynum = 0;
		for (int i = 0; i < 3; ++i) {
			if (xppac_xflag & (1 << i)) ++channel.ppac_xnum;
			if (xppac_yflag & (1 << i)) ++channel.ppac_ynum;
		}
		if (channel.ppac_xnum < 1 || channel.ppac_ynum < 1) {
			channel.valid |= 4;
			continue;
		}
		if (channel.ppac_xnum >= 2 && channel.ppac_ynum >= 2) {
			channel.ppac_valid = true;
		} else {
			channel.ppac_valid = false;
			channel.valid |= 4;
		}

		// get target position
		// correct PPAC
		for (int i = 0; i < 3; ++i) {
			channel.ppac_x[i] = xppac_event.x[i] - using_ppac_correct[0][i];
			channel.ppac_y[i] = xppac_event.y[i] - using_ppac_correct[1][i];
		}
		// track
		double xk, yk;
		TrackMultiplePpac(
			channel.ppac_xflag, using_ppac_xz, channel.ppac_x, xk, channel.tx
		);
		TrackMultiplePpac(
			channel.ppac_yflag, using_ppac_yz, channel.ppac_y, yk, channel.ty
		);

		// get recoil direction
		ROOT::Math::XYZVector recoil_direction(
			channel.recoil_x - channel.tx,
			channel.recoil_y - channel.ty,
			channel.recoil_z
		);
		recoil_direction = recoil_direction.Unit();
		// get recoil energy
		if (channel.tafcsi_valid) {
			channel.recoil_kinetic = taf_events[channel.taf_index].energy[0];
		} else {
			channel.recoil_kinetic = h2_calculator.Energy(
				recoil_direction.Theta(),
				taf_events[channel.taf_index].energy[0],
				tafd_thickness[channel.taf_index]
			);
		}
		channel.recoil_energy = channel.recoil_kinetic + recoil_mass;
		channel.recoil_momentum = MomentumFromKinetic(
			recoil_mass, channel.recoil_kinetic
		);
		// get recoil position
		channel.recoil_x = taf_events[channel.taf_index].x[0];
		channel.recoil_y = taf_events[channel.taf_index].y[0];
		channel.recoil_z = taf_events[channel.taf_index].z[0];
		ROOT::Math::XYZVector recoil_p = recoil_direction * channel.recoil_momentum;

		// get fragment energy
		for (int i = 0; i < 2; ++i) {
			channel.fragment_kinetic[i] = t0_event.energy[frag_index[i]];
			channel.fragment_energy[i] =
				channel.fragment_kinetic[i] + fragment_mass[i];
			channel.fragment_momentum[i] = MomentumFromKinetic(
				fragment_mass[i], channel.fragment_kinetic[i]
			);
		}
		// get fragment position
		for (int i = 0; i < 2; ++i) {
			channel.fragment_x[i] = t0_event.x[frag_index[i]];
			channel.fragment_y[i] = t0_event.y[frag_index[i]];
			channel.fragment_z[i] = t0_event.z[frag_index[i]];
		}
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

		// get parent momentum
		ROOT::Math::XYZVector parent_p = fragment_p[0] + fragment_p[1];
		channel.parent_momentum = parent_p.R();

		// get beam momentum and energy
		ROOT::Math::XYZVector beam_p =
			fragment_p[0] + fragment_p[1] + recoil_p;
		channel.beam_momentum = beam_p.R();
		channel.beam_kinetic = KineticFromMomentum(
			beam_mass, channel.beam_momentum
		);
		channel.beam_energy =
			pow(beam_mass, 2.0) + pow(channel.beam_momentum, 2.0);

		// other information
		channel.run = run;
		channel.entry = entry;

		// extra TAF information
		channel.tafd_front_strip =
			taf_tele_events[channel.taf_index].front_strip[0];
		channel.tafd_energy =
			taf_tele_events[channel.taf_index].energy[0][0];
		channel.t0_layer[0] = t0_event.layer[frag_index[0]];
		channel.t0_layer[1] = t0_event.layer[frag_index[1]];

		opt.Fill();
	}
	// show finish
	if (simulation) {
		printf("\b\b\b\b100%%\n");
	}
}

void PrintUsage(const char* name) {
	std::cout << "Usage: " << name << " [options] [recoil]\n"
		<< "  recoil        Recoil particle mass, default is 2(2H).\n"
		<< "Options:\n"
		<< "  -s            Process simulated data.\n"
		<< "  -c            Use C target.\n";
}

int main(int argc, char **argv) {
	if (argc > 2) {
		PrintUsage(argv[0]);
		return -1;
	}
	bool simulation = false;
	int recoil_mass_number = 2;
	std::string target = "cd2";
	if (argc == 2) {
		if (argv[1][0] == '-' && argv[1][1] == 'h') {
			PrintUsage(argv[0]);
			return 0;
		} else if (argv[1][0] == '-' && argv[1][1] == 's') {
			simulation = true;
		} else if (argv[1][0] == '-' && argv[1][1] == 'c') {
			target = "c";
		} else {
			recoil_mass_number = atoi(argv[1]);
		}
	}
	if (recoil_mass_number != 1 && recoil_mass_number != 2) {
		std::cerr << "Error: Invalid recoil mass number, "
			<< recoil_mass_number << "\n";
		return -1;
	}
	if (recoil_mass_number == 1) recoil_mass = mass_1h;

	if (simulation) std::cout << "Use simulated data.\n";
	std::cout << "Recoil particle is " << recoil_mass_number << "H\n";

	// output file name
	TString output_file_name = TString::Format(
		"%s%sC14-10Be-4He-%dH-v2%s%s.root",
		kGenerateDataPath,
		kChannelDir,
		recoil_mass_number,
		target == "c" ? "-c" : "",
		simulation ? "-sim" : ""
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "channel");
	// output event
	ChannelV2Event channel;
	// setup output branches
	channel.SetupOutput(&opt);

	// initialize
	channel.beam_charge = 6;
	channel.beam_mass = 14;
	channel.parent_charge = 6;
	channel.parent_mass = 14;
	channel.recoil_charge = 1;
	channel.recoil_mass = recoil_mass_number;
	channel.fragment_num = 2;
	channel.fragment_charge[0] = 4;
	channel.fragment_mass[0] = 10;
	channel.fragment_charge[1] = 2;
	channel.fragment_mass[1] = 4;

	// statistics
	int taf_conflict_number = 0;

	// calculator
	elc::CsiEnergyCalculator h2_calculator("2H");

	if (simulation) {
		CoincideSingleRun(2, true, opt, channel, h2_calculator, taf_conflict_number);
	} else if (target == "cd2") {
		printf("Processing run 000");
		for (int run = 618; run <= 716; ++run) {
			if (run == 628) continue;
			if (run > 652 && run < 675) continue;
			if (run > 716 && run < 739) continue;

			printf("\b\b\b%3d", run);
			fflush(stdout);

			CoincideSingleRun(run, false, opt, channel, h2_calculator, taf_conflict_number);
		}
		// show finish
		printf("\n");
	} else if (target == "c") {
		printf("Processing run 000");
		for (int run = 656; run <= 732; ++run) {
			if (run == 723) continue;
			if (run > 657 && run < 668) continue;
			if (run > 669 && run < 717) continue;

			printf("\b\b\b%3d", run);
			fflush(stdout);

			CoincideSingleRun(run, false, opt, channel, h2_calculator, taf_conflict_number);
		}
		// show finish
		printf("\n");
	}

	std::cout << "Conflict TAF events " << taf_conflict_number << "\n";

	// save
	opf.cd();
	opt.Write();
	opf.Close();

	return 0;
}