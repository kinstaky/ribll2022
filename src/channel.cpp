#include "include/channel.h"

#include <iostream>
#include <map>
#include <exception>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/particle_event.h"
#include "include/event/channel_event.h"

namespace ribll {

const std::map<std::string, unsigned short> element{
	{"H", 1},
	{"He", 2},
	{"Li", 3},
	{"Be", 4},
	{"B", 5},
	{"C", 6}
};

double IonMass(unsigned short charge, unsigned short mass) {
	constexpr double electron_mass = 5.48579909065e-4;
	if (charge == 1) {
		if (mass == 1) return 1.0072764520;
		else if (mass == 2) return 2.0135531980;
		else if (mass == 3) return 3.0155007014;
	} else if (charge == 2) {
		// He
		if (mass == 3) return 3.0149321622;
		else if (mass == 4) return 4.0015060943;
		else if (mass == 6) return 6.0177887354;
	} else if (charge == 3) {
		// Li
		if (mass == 6) return 6.0134771477;
		else if (mass == 7) return  7.0143576949;
	} else if (charge == 4) {
		// Be
		if (mass == 7) return 7.0147343973;
		else if (mass == 8) return 8.0031108;
		else if (mass == 9) return 9.0099887420;
		else if (mass == 10) return 10.0113403769;
	} else if (charge == 5) {
		// B
		if (mass == 10) return 10.0101939628;
		else if (mass == 11) return 11.0065622673;
	} else if (charge == 6) {
		// C
		if (mass == 12) return 11.9967085206;
		else if (mass == 13) return 13.0000633559;
		else if (mass == 14) return 13.9999505089;
		else if (mass == 15) return 15.0073077289;
	}
	return double(mass) - charge * electron_mass;
}


/// @brief calculate momentum from energy considering relative effects
/// @param[in] kinetic_energy kinetic energy of particle, in MeV
/// @param[in] mass mass number of particle
/// @returns momentum value of particle
///
double MomentumFromEnergy(
	double kinetic_energy,
	unsigned short charge,
	unsigned short mass
) {
	// atomic mass constant
	constexpr double u = 931.494;
	double ion_mass = IonMass(charge, mass) * u;
	// double energy = kinetic_energy + mass * u;
	// double p = sqrt(energy * energy - mass*mass*u*u);
	return sqrt(kinetic_energy * (kinetic_energy + 2.0*ion_mass));
}


double EnergyFromMomentum(
	double momentum,
	unsigned short charge,
	unsigned short mass
) {
	// atomic mass constant
	constexpr double u = 931.494;
	// accurate mass in MeV/c^2
	double ion_mass = IonMass(charge, mass) * u;

	return sqrt(momentum*momentum + ion_mass*ion_mass) - ion_mass;
}


//-----------------------------------------------------------------------------
// 								Channel
//-----------------------------------------------------------------------------

Channel::Channel(unsigned int run)
: run_(run) {
}


int Channel::Coincide() {
	std::cerr << "Error: Channel::Coincide not implemented yet.\n";
	return -1;
}

//-----------------------------------------------------------------------------
// 							Be8ToTwoAlphaChannel
//-----------------------------------------------------------------------------

Be8ToTwoAlphaChannel::Be8ToTwoAlphaChannel(unsigned int run)
: Channel(run) {
}


int Be8ToTwoAlphaChannel::Coincide() {
	// t0 particle file name
	TString t0_file_name;
	t0_file_name.Form(
		"%s%st0-particle-%04d.root",
		kGenerateDataPath,
		kParticleDir,
		run_
	);
	// t0 particle file
	TFile t0_file(t0_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)t0_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< t0_file_name << " failed.\n";
		return -1;
	}
	// add XIA PPAC friend
	ipt->AddFriend("xppac=tree", TString::Format(
		"%s%sxppac-particle-%04d.root",
		kGenerateDataPath,
		kParticleDir,
		run_
	));
	// input T0 particle events
	ParticleEvent t0;
	// PPAC particle event
	ParticleEvent xppac;
	//setup input branches
	t0.SetupInput(ipt);
	xppac.SetupInput(ipt, "xppac.");

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%sBe8-4He-4He-%04u.root",
		kGenerateDataPath,
		kChannelDir,
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "channel");
	// output channel event
	ChannelEvent channel;
	// setup output branches
	channel.SetupOutput(&opt);

	// total valid events
	long long total = 0;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Coinciding particles   0%%");
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

		int valid_alpha = 0;
		unsigned short alpha_index[16];
		for (unsigned short i = 0; i < t0.num; ++i) {
			if (t0.charge[i] == 2 && t0.mass[i] == 4) {
				alpha_index[valid_alpha] = i;
				++valid_alpha;
			}
		}
		if (valid_alpha != 2) continue;
		// check PPAC tracking
		if (xppac.num != 4) continue;
		if (pow(xppac.x[3]+3.3, 2.0)+pow(xppac.y[3]-1.0, 2.0) > 225.0) {
			continue;
		}

		// fill channel particle number
		channel.num = 2;
		// momentums
		std::vector<ROOT::Math::XYZVector> p;
		// fill particles from T0
		for (unsigned short i = 0; i < channel.num; ++i) {
			channel.daughter_charge[i] = t0.charge[alpha_index[i]];
			channel.daughter_mass[i] = t0.mass[alpha_index[i]];
			channel.daughter_energy[i] = t0.energy[alpha_index[i]];
			channel.daughter_time[i] = t0.time[alpha_index[i]];
			// calculate momentum
			double momentum = MomentumFromEnergy(
				channel.daughter_energy[i],
				channel.daughter_charge[i],
				channel.daughter_mass[i]
			);
			p.emplace_back(
				t0.x[alpha_index[i]] - xppac.x[3],
				t0.y[alpha_index[i]] - xppac.y[3],
				t0.z[alpha_index[i]] - xppac.z[3]
			);
			p[i] = p[i].Unit() * momentum;
			channel.daughter_px[i] = p[i].X();
			channel.daughter_py[i] = p[i].Y();
			channel.daughter_pz[i] = p[i].Z();
			channel.daughter_r[i] = p[i].R();
			channel.daughter_theta[i] = p[i].Theta();
			channel.daughter_phi[i] = p[i].Phi();
		}

		// fill beam from PPAC
		channel.beam_energy = 375.0;
		channel.beam_time = xppac.time[3];
		// beam p
		ROOT::Math::XYZVector bp(xppac.px[3], xppac.py[3], xppac.pz[3]);
		double bp_value = MomentumFromEnergy(channel.beam_energy, 6, 14);
		bp *= bp_value;
		channel.beam_px = bp.X();
		channel.beam_py = bp.Y();
		channel.beam_pz = bp.Z();
		channel.beam_r = bp.R();
		channel.beam_theta = bp.Theta();
		channel.beam_phi = bp.Phi();

		// rebuild recoil particle
		channel.recoil = 0;

		// fill parent nuclear
		channel.parent_charge = 4;
		channel.parent_mass = 8;
		// parent p
		ROOT::Math::XYZVector pp(0.0, 0.0, 0.0);
		for (unsigned short i = 0; i < channel.num; ++i) {
			pp += p[i];
		}
		channel.parent_energy = 0.0;
		for (unsigned short i = 0; i < channel.num; ++i) {
			channel.parent_energy += channel.daughter_energy[i]
				+ IonMass(
					channel.daughter_charge[i],
					channel.daughter_mass[i]
				) * 931.494;
		}
		// parent mass
		double pm = sqrt(channel.parent_energy * channel.parent_energy - pp.Mag2());
		// parent kinematic energy
		channel.parent_energy -= pm;

		// fill parent momentum
		channel.parent_px = pp.X();
		channel.parent_py = pp.Y();
		channel.parent_pz = pp.Z();
		channel.parent_r = pp.R();
		channel.parent_theta = pp.Theta();
		channel.parent_phi = pp.Phi();

		// entry
		channel.entry = entry;
		// taf index
		channel.taf_index = -1;

		++total;
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	std::cout << "Total valid event is " << total << "\n";

	// save tree
	opt.Write();
	// close file
	t0_file.Close();
	opf.Close();
	return 0;
}



//-----------------------------------------------------------------------------
// 							C12ToThreeAlphaChannel
//-----------------------------------------------------------------------------

C12ToThreeAlphaChannel::C12ToThreeAlphaChannel(unsigned int run)
: Channel(run) {
}


int C12ToThreeAlphaChannel::Coincide() {
	// t0 particle file name
	TString t0_file_name;
	t0_file_name.Form(
		"%s%st0-particle-%04d.root",
		kGenerateDataPath,
		kParticleDir,
		run_
	);
	// t0 particle file
	TFile t0_file(t0_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)t0_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< t0_file_name << " failed.\n";
		return -1;
	}
	// add XIA PPAC friend
	ipt->AddFriend("xppac=tree", TString::Format(
		"%s%sxppac-particle-%04d.root",
		kGenerateDataPath,
		kParticleDir,
		run_
	));
	// input T0 particle events
	ParticleEvent t0;
	// PPAC particle event
	ParticleEvent xppac;
	//setup input branches
	t0.SetupInput(ipt);
	xppac.SetupInput(ipt, "xppac.");

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%sC12-4He-4He-4He-%04u.root",
		kGenerateDataPath,
		kChannelDir,
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "channel");
	// output channel event
	ChannelEvent channel;
	// setup output branches
	channel.SetupOutput(&opt);

	// total valid events
	long long total = 0;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Coinciding particles   0%%");
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

		int valid_alpha = 0;
		unsigned short alpha_index[16];
		for (unsigned short i = 0; i < t0.num; ++i) {
			if (t0.charge[i] == 2 && t0.mass[i] == 4) {
				alpha_index[valid_alpha] = i;
				++valid_alpha;
			}
		}
		if (valid_alpha != 3) continue;
		// check PPAC tracking
		if (xppac.num != 4) continue;
		if (pow(xppac.x[3]+3.3, 2.0)+pow(xppac.y[3]-1.0, 2.0) > 225.0) {
			continue;
		}

		// fill channel particle number
		channel.num = 3;
		// momentums
		std::vector<ROOT::Math::XYZVector> p;
		// fill particles from T0
		for (unsigned short i = 0; i < channel.num; ++i) {
			channel.daughter_charge[i] = t0.charge[alpha_index[i]];
			channel.daughter_mass[i] = t0.mass[alpha_index[i]];
			channel.daughter_energy[i] = t0.energy[alpha_index[i]];
			channel.daughter_time[i] = t0.time[alpha_index[i]];
			// calculate momentum
			double momentum = MomentumFromEnergy(
				channel.daughter_energy[i],
				channel.daughter_charge[i],
				channel.daughter_mass[i]
			);
			p.emplace_back(
				t0.x[alpha_index[i]] - xppac.x[3],
				t0.y[alpha_index[i]] - xppac.y[3],
				t0.z[alpha_index[i]] - xppac.z[3]
			);
			p[i] = p[i].Unit() * momentum;
			channel.daughter_px[i] = p[i].X();
			channel.daughter_py[i] = p[i].Y();
			channel.daughter_pz[i] = p[i].Z();
			channel.daughter_r[i] = p[i].R();
			channel.daughter_theta[i] = p[i].Theta();
			channel.daughter_phi[i] = p[i].Phi();
		}

		// fill beam from PPAC
		channel.beam_energy = 375.0;
		channel.beam_time = xppac.time[3];
		// beam p
		ROOT::Math::XYZVector bp(xppac.px[3], xppac.py[3], xppac.pz[3]);
		double bp_value = MomentumFromEnergy(channel.beam_energy, 6, 14);
		bp *= bp_value;
		channel.beam_px = bp.X();
		channel.beam_py = bp.Y();
		channel.beam_pz = bp.Z();
		channel.beam_r = bp.R();
		channel.beam_theta = bp.Theta();
		channel.beam_phi = bp.Phi();

		// rebuild recoil particle
		channel.recoil = 0;

		// fill parent nuclear
		channel.parent_charge = 6;
		channel.parent_mass = 12;
		// parent p
		ROOT::Math::XYZVector pp(0.0, 0.0, 0.0);
		for (unsigned short i = 0; i < channel.num; ++i) {
			pp += p[i];
		}
		channel.parent_energy = 0.0;
		for (unsigned short i = 0; i < channel.num; ++i) {
			channel.parent_energy += channel.daughter_energy[i]
				+ IonMass(
					channel.daughter_charge[i],
					channel.daughter_mass[i]
				) * 931.494;
		}
		// parent mass
		double pm = sqrt(channel.parent_energy * channel.parent_energy - pp.Mag2());
		// parent kinematic energy
		channel.parent_energy -= pm;

		// fill parent momentum
		channel.parent_px = pp.X();
		channel.parent_py = pp.Y();
		channel.parent_pz = pp.Z();
		channel.parent_r = pp.R();
		channel.parent_theta = pp.Theta();
		channel.parent_phi = pp.Phi();

		// entry
		channel.entry = entry;
		// taf index
		channel.taf_index = -1;

		++total;
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	std::cout << "Total valid event is " << total << "\n";

	// save tree
	opt.Write();
	// close file
	t0_file.Close();
	opf.Close();
	return 0;
}


//-----------------------------------------------------------------------------
// 							C14ToBe10He4TwoBodyChannel
//-----------------------------------------------------------------------------

C14ToBe10He4TwoBodyChannel::C14ToBe10He4TwoBodyChannel(unsigned int run)
: Channel(run) {
}


int C14ToBe10He4TwoBodyChannel::Coincide() {
	// t0 particle file name
	TString t0_file_name;
	t0_file_name.Form(
		"%s%st0-particle-%04d.root",
		kGenerateDataPath,
		kParticleDir,
		run_
	);
	// t0 particle file
	TFile t0_file(t0_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)t0_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< t0_file_name << " failed.\n";
		return -1;
	}
	// add XIA PPAC friend
	ipt->AddFriend("xppac=tree", TString::Format(
		"%s%sxppac-particle-%04d.root",
		kGenerateDataPath,
		kParticleDir,
		run_
	));
	// input T0 particle events
	ParticleEvent t0;
	// PPAC particle event
	ParticleEvent xppac;
	//setup input branches
	t0.SetupInput(ipt);
	xppac.SetupInput(ipt, "xppac.");

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%sC14-10Be-4He-%04u.root",
		kGenerateDataPath,
		kChannelDir,
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "channel");
	// output channel event
	ChannelEvent channel;
	// setup output branches
	channel.SetupOutput(&opt);

	// total valid events
	long long total = 0;
	long long conflict_t0 = 0;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Coinciding particles   0%%");
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

		// check T0 event
		int t0_status = 0;
		size_t index[4];
		bool t0_valid = true;
		for (unsigned short i = 0; i < t0.num; ++i) {
			if (t0.charge[i] == 4 && t0.mass[i] == 10) {
				// find 10Be
				if ((t0_status & 0x1) == 0) {
					index[0] = i;
					t0_status |= 0x1;
				} else {
					t0_valid = false;
				}
			} else if (t0.charge[i] == 2 && t0.mass[i] == 4) {
				// find 4He
				if ((t0_status & 0x2) == 0) {
					index[1] = i;
					t0_status |= 0x2;
				} else {
					t0_valid = false;
				}
			}
		}
		if (!t0_valid && t0_status == 0x3) ++conflict_t0;
		if (!t0_valid || t0_status != 0x3) continue;

		// check PPAC tracking
		if (xppac.num != 4) continue;
		if (pow(xppac.x[3]+3.3, 2.0)+pow(xppac.y[3]-1.0, 2.0) > 225.0) {
			continue;
		}

		// fill channel particle number
		channel.num = 2;
		// momentums
		std::vector<ROOT::Math::XYZVector> p;
		// fill particles from T0
		for (unsigned short i = 0; i < channel.num; ++i) {
			channel.daughter_charge[i] = t0.charge[index[i]];
			channel.daughter_mass[i] = t0.mass[index[i]];
			channel.daughter_energy[i] = t0.energy[index[i]];
			channel.daughter_time[i] = t0.time[index[i]];
			// calculate momentum
			double momentum = MomentumFromEnergy(
				channel.daughter_energy[i],
				channel.daughter_charge[i],
				channel.daughter_mass[i]
			);
			p.emplace_back(
				t0.x[index[i]] - xppac.x[3],
				t0.y[index[i]] - xppac.y[3],
				t0.z[index[i]] - xppac.z[3]
			);
			p[i] = p[i].Unit() * momentum;
			channel.daughter_px[i] = p[i].X();
			channel.daughter_py[i] = p[i].Y();
			channel.daughter_pz[i] = p[i].Z();
			channel.daughter_r[i] = p[i].R();
			channel.daughter_theta[i] = p[i].Theta();
			channel.daughter_phi[i] = p[i].Phi();
		}

		// fill beam from PPAC
		channel.beam_energy = 375.0;
		channel.beam_time = xppac.time[3];
		// beam p
		ROOT::Math::XYZVector bp(xppac.px[3], xppac.py[3], xppac.pz[3]);
		double bp_value = MomentumFromEnergy(channel.beam_energy, 6, 14);
		bp *= bp_value;
		channel.beam_px = bp.X();
		channel.beam_py = bp.Y();
		channel.beam_pz = bp.Z();
		channel.beam_r = bp.R();
		channel.beam_theta = bp.Theta();
		channel.beam_phi = bp.Phi();

		// rebuild recoil particle
		channel.recoil = 0;
		ROOT::Math::XYZVector rp = bp - p[0] - p[1];
		channel.recoil_charge = 1;
		channel.recoil_mass = 2;
		channel.recoil_energy = EnergyFromMomentum(
			rp.R(), channel.recoil_charge, channel.recoil_mass
		);
		channel.recoil_px = rp.X();
		channel.recoil_py = rp.Y();
		channel.recoil_pz = rp.Z();
		channel.recoil_r = rp.R();
		channel.recoil_theta = rp.Theta();
		channel.recoil_phi = rp.Phi();

		// fill parent nuclear
		channel.parent_charge = 6;
		channel.parent_mass = 14;
		// parent p
		ROOT::Math::XYZVector pp = bp;
		channel.parent_energy = channel.beam_energy;
		// fill parent momentum
		channel.parent_px = pp.X();
		channel.parent_py = pp.Y();
		channel.parent_pz = pp.Z();
		channel.parent_r = pp.R();
		channel.parent_theta = pp.Theta();
		channel.parent_phi = pp.Phi();

		// entry
		channel.entry = entry;
		// taf index
		channel.taf_index = -1;

		++total;
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	std::cout << "Total valid event is " << total << "\n"
		<< "Conflict T0 event " << conflict_t0 << "\n";

	// save tree
	opt.Write();
	// close file
	t0_file.Close();
	opf.Close();
	return 0;
}


//-----------------------------------------------------------------------------
// 							C14ToBe10He4ThreeBodyChannel
//-----------------------------------------------------------------------------

C14ToBe10He4ThreeBodyChannel::C14ToBe10He4ThreeBodyChannel(unsigned int run)
: Channel(run) {
}


int C14ToBe10He4ThreeBodyChannel::Coincide() {
	// t0 particle file name
	TString t0_file_name;
	t0_file_name.Form(
		"%s%st0-particle-ta-%04d.root",
		kGenerateDataPath,
		kParticleDir,
		run_
	);
	// t0 particle file
	TFile t0_file(t0_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)t0_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< t0_file_name << " failed.\n";
		return -1;
	}
	// add TAF friends
	for (int i = 0 ; i < 6; ++i) {
		ipt->AddFriend(
			TString::Format("taf%d=tree", i),
			TString::Format(
				"%s%staf%d-particle-ta-%04d.root",
				kGenerateDataPath,
				kParticleDir,
				i,
				run_
			)
		);
	}
	// add XIA PPAC friend
	ipt->AddFriend("xppac=tree", TString::Format(
		"%s%sxppac-particle-ta-%04d.root",
		kGenerateDataPath,
		kParticleDir,
		run_
	));
	// input T0 particle events
	ParticleEvent t0;
	// input TAF particle events
	ParticleEvent taf[6];
	// input XIA PPAC particle events
	ParticleEvent xppac;
	// setup input branches
	t0.SetupInput(ipt);
	for (int i = 0; i < 6; ++i) {
		taf[i].SetupInput(ipt, "taf"+std::to_string(i)+".");
	}
	xppac.SetupInput(ipt, "xppac.");


	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%sC14-10Be-4He-2H-%04u.root",
		kGenerateDataPath,
		kChannelDir,
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "channel");
	// output channel event
	ChannelEvent channel;
	// setup output branches
	channel.SetupOutput(&opt);

	// total valid events
	long long total = 0;
	long long conflict_taf = 0;
	long long conflict_t0 = 0;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Coinciding particles   0%%");
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

		// check T0 event
		int t0_status = 0;
		size_t index[4];
		bool t0_valid = true;
		for (unsigned short i = 0; i < t0.num; ++i) {
			if (t0.charge[i] == 4 && t0.mass[i] == 10) {
				// find 10Be
				if ((t0_status & 0x1) == 0) {
					index[0] = i;
					t0_status |= 0x1;
				} else {
					t0_valid = false;
				}
			} else if (t0.charge[i] == 2 && t0.mass[i] == 4) {
				// find 4He
				if ((t0_status & 0x2) == 0) {
					index[1] = i;
					t0_status |= 0x2;
				} else {
					t0_valid = false;
				}
			}
		}
		if (!t0_valid && t0_status == 0x3) ++conflict_t0;
		if (!t0_valid || t0_status != 0x3) continue;
		// check TAF
		int taf_index = -1;
		int valid = 0;
		for (int i = 0; i < 6; ++i) {
			if (
				taf[i].num == 1
				&& taf[i].charge[0] == 1
				&& taf[i].mass[0] == 2
				&& taf[i].energy[0] > -9e4
			) {
				++valid;
				taf_index = i;
			}
		}
		// jump events without 2H
		if (valid == 0) continue;
		// jump events with more than one 2H
		if (valid > 1) {
			++conflict_taf;
			continue;
		}
		// check PPAC tracking
		if (xppac.num != 4) continue;
		if (pow(xppac.x[3]+3.3, 2.0)+pow(xppac.y[3]-1.0, 2.0) > 225.0) {
			continue;
		}

		// fill channel particle number
		channel.num = 2;
		// momentums
		std::vector<ROOT::Math::XYZVector> p;
		// fill particles from T0
		for (unsigned short i = 0; i < channel.num; ++i) {
			channel.daughter_charge[i] = t0.charge[index[i]];
			channel.daughter_mass[i] = t0.mass[index[i]];
			channel.daughter_energy[i] = t0.energy[index[i]];
			channel.daughter_time[i] = t0.time[index[i]];
			// calculate momentum
			double momentum = MomentumFromEnergy(
				channel.daughter_energy[i],
				channel.daughter_charge[i],
				channel.daughter_mass[i]
			);
			p.emplace_back(
				t0.x[index[i]] - xppac.x[3],
				t0.y[index[i]] - xppac.y[3],
				t0.z[index[i]] - xppac.z[3]
			);
			p[i] = p[i].Unit() * momentum;
			channel.daughter_px[i] = p[i].X();
			channel.daughter_py[i] = p[i].Y();
			channel.daughter_pz[i] = p[i].Z();
			channel.daughter_r[i] = p[i].R();
			channel.daughter_theta[i] = p[i].Theta();
			channel.daughter_phi[i] = p[i].Phi();
		}

		// fill particles from TAF
		channel.recoil = 1;
		channel.recoil_charge = taf[taf_index].charge[0];
		channel.recoil_mass = taf[taf_index].mass[0];
		channel.recoil_energy = taf[taf_index].energy[0];
		channel.recoil_time = taf[taf_index].time[0];
		// calculate momentum
		double momentum = MomentumFromEnergy(
			channel.recoil_energy,
			channel.recoil_charge,
			channel.recoil_mass
		);
		// recoil p
		ROOT::Math::XYZVector rp(
			taf[taf_index].x[0] - xppac.x[3],
			taf[taf_index].y[0] - xppac.y[3],
			taf[taf_index].z[0] - xppac.z[3]
		);
		rp = rp.Unit() * momentum;
		channel.recoil_px = rp.X();
		channel.recoil_py = rp.Y();
		channel.recoil_pz = rp.Z();
		channel.recoil_r = rp.R();
		channel.recoil_theta = rp.Theta();
		channel.recoil_phi = rp.Phi();

		// fill beam from PPAC
		channel.beam_energy = 375.0;
		channel.beam_time = xppac.time[3];
		ROOT::Math::XYZVector bp(xppac.px[3], xppac.py[3], xppac.pz[3]);
		double bp_value = MomentumFromEnergy(channel.beam_energy, 6, 14);
		bp = bp.Unit() * bp_value;
		channel.beam_px = bp.X();
		channel.beam_py = bp.Y();
		channel.beam_pz = bp.Z();
		channel.beam_r = bp.R();
		channel.beam_theta = bp.Theta();
		channel.beam_phi = bp.Phi();

		// fill parent nuclear
		channel.parent_charge = 6;
		channel.parent_mass = 14;
		// parent p
		ROOT::Math::XYZVector pp(0.0, 0.0, 0.0);
		for (unsigned short i = 0; i < channel.num; ++i) {
			pp += p[i];
		}
		pp += rp;
		channel.parent_energy = EnergyFromMomentum(
			pp.R(), channel.parent_charge, channel.parent_mass
		);
		// for (unsigned short i = 0; i < channel.num; ++i) {
		// 	channel.parent_energy += channel.daughter_energy[i]
		// 		+ IonMass(channel.daughter_charge[i], channel.daughter_mass[i]) * 931.494;
		// }
		// double pm = sqrt(channel.parent_energy * channel.parent_energy - pp.Mag2());
		// channel.parent_energy -= pm;

		channel.parent_px = pp.X();
		channel.parent_py = pp.Y();
		channel.parent_pz = pp.Z();
		channel.parent_r = pp.R();
		channel.parent_theta = pp.Theta();
		channel.parent_phi = pp.Phi();

		// entry
		channel.entry = entry;
		// taf index
		channel.taf_index = taf_index;

		++total;
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	std::cout << "Total valid event is " << total << "\n"
		<< "Conflict TAF event " << conflict_taf << "\n"
		<< "Conflcit T0 event " << conflict_t0 << "\n";

	// save tree
	opt.Write();
	// close file
	t0_file.Close();
	opf.Close();
	return 0;
}


}	// namespace ribll