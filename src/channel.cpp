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
			// status
			channel.status[i] = t0.status[alpha_index[i]];
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
		// parent kinetic energy
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
			// other information
			channel.status[i] = t0.status[alpha_index[i]];
		}

		// fill beam from PPAC
		channel.beam_energy = 375.0;
		channel.beam_time = xppac.time[3];
		// beam p
		ROOT::Math::XYZVector bp(xppac.px[3], xppac.py[3], xppac.pz[3]);
		double bp_value = MomentumFromEnergy(channel.beam_energy, 6, 12);
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
		for (int i = 0; i < channel.num; ++i) {
			channel.parent_energy += channel.daughter_energy[i]
				+ IonMass(
					channel.daughter_charge[i],
					channel.daughter_mass[i]
				) * 931.494;
		}
		// parent mass
		double pm = sqrt(channel.parent_energy * channel.parent_energy - pp.Mag2());
		// parent kinetic energy
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
		channel.beam_energy = 392.0;
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

C14ToBe10He4ThreeBodyChannel::C14ToBe10He4ThreeBodyChannel(
	unsigned int run,
	unsigned short recoil_mass,
	bool simulate
)
: Channel(run)
, recoil_mass_(recoil_mass)
, simulate_(simulate) {
}


int C14ToBe10He4ThreeBodyChannel::Coincide() {
	const int be_mass_number = 10;
	// t0 particle file name
	TString t0_file_name;
	t0_file_name.Form(
		"%s%st0-particle-%sta-%04d.root",
		kGenerateDataPath,
		kParticleDir,
		simulate_ ? "sim-" : "",
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
				"%s%staf%d-particle-%sta-v2-%04d.root",
				kGenerateDataPath,
				kParticleDir,
				i,
				simulate_ ? "sim-" : "",
				run_
			)
		);
	}
	// add XIA PPAC friend
	ipt->AddFriend("xppac=tree", TString::Format(
		"%s%sxppac-particle-%sta-%04d.root",
		kGenerateDataPath,
		kParticleDir,
		simulate_ ? "sim-" : "",
		run_
	));
	// add VME PPAC friend
	ipt->AddFriend("vppac=tree", TString::Format(
		"%s%svppac-particle-%sta-%04d.root",
		kGenerateDataPath,
		kParticleDir,
		simulate_ ? "sim-" : "",
		run_
	));
	// input T0 particle events
	ParticleEvent t0;
	bool hole[4];
	// input TAF particle events
	ParticleEvent taf[6];
	// input XIA PPAC particle events
	ParticleEvent xppac;
	unsigned short xppac_xflag;
	unsigned short xppac_yflag;
	// input VME PPAC particle event
	ParticleEvent vppac;
	unsigned short vppac_xflag;
	unsigned short vppac_yflag;
	// setup input branches
	t0.SetupInput(ipt);
	ipt->SetBranchAddress("hole", hole);
	for (int i = 0; i < 6; ++i) {
		taf[i].SetupInput(ipt, "taf"+std::to_string(i)+".");
	}
	xppac.SetupInput(ipt, "xppac.");
	ipt->SetBranchAddress("xppac.xflag", &xppac_xflag);
	ipt->SetBranchAddress("xppac.yflag", &xppac_yflag);
	vppac.SetupInput(ipt, "vppac.");
	ipt->SetBranchAddress("vppac.xflag", &vppac_xflag);
	ipt->SetBranchAddress("vppac.yflag", &vppac_yflag);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%sC14-%dBe-4He-%dH-%s%04u.root",
		kGenerateDataPath,
		kChannelDir,
		be_mass_number,
		recoil_mass_,
		simulate_ ? "sim-" : "",
		run_
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "channel");
	// output channel event
	ChannelEvent channel;
	int t0_index[8];
	double q_value;
	int csi_index;
	int ppac_flag;
	// PPAC xhit and yhit
	int xhit, yhit;
	// XPPAC num
	int xnum;
	// setup output branches
	channel.SetupOutput(&opt);
	opt.Branch("t0_index", t0_index, "t0_index[num]/I");
	opt.Branch("q", &q_value, "q/D");
	opt.Branch("csi_index", &csi_index, "ci/I");
	opt.Branch("ppac_flag", &ppac_flag, "pflag/I");
	opt.Branch("ppac_xhit", &xhit, "pxhit/I");
	opt.Branch("ppac_yhit", &yhit, "pyhit/I");
	opt.Branch("xppac_num", &xnum, "xnum/I");

	// total valid events
	long long total = 0;
	long long conflict_taf = 0;
	long long conflict_t0 = 0;
	long long possible_2h = 0;

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

		// check T0 and find 10Be + 4He
		int t0_status = 0;
		size_t index[4];
		bool t0_valid = true;
		for (unsigned short i = 0; i < t0.num; ++i) {
			// ignore bad strips
			// if (hole[i]) continue;
			if (t0.charge[i] == 4 && t0.mass[i] == be_mass_number) {
				// find 10Be
				if ((t0_status & 0x1) == 0) {
					index[0] = i;
					t0_status |= 0x1;
					t0.mass[i] = be_mass_number;
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
		// check TAF and find 2H
		channel.taf_index = -1;
		int valid = 0;
		// loop TAFs and search 2H with dE-E PID
		for (int i = 0; i < 6; ++i) {
			if (
				taf[i].num == 1
				&& taf[i].charge[0] == 1
				&& taf[i].mass[0] == recoil_mass_
				&& taf[i].energy[0] > -9e4
			) {
				++valid;
				channel.taf_index = i;
			}
		}
		// if (valid == 0) {
		// 	// loop TAFs again if not found, try to find 2H stopped in TAFD
		// 	for (int i = 0; i < 6; ++i) {
		// 		// charge==201 means stopped in ADSSD
		// 		if (
		// 			taf[i].num == 1
		// 			&& taf[i].charge[0] == 201
		// 			&& taf[i].energy[0] > -9e4
		// 		) {
		// 			++valid;
		// 		}
		// 	}
		// 	if (valid > 0) ++possible_2h;
		// }

		// check valid
		if (!t0_valid && t0_status == 0x3) ++conflict_t0;
		if (!t0_valid || t0_status != 0x3) {
			if (simulate_) {
				if (valid == 0) {
					channel.taf_index = -6;
				} else {
					channel.taf_index = -2;
				}
				opt.Fill();
			}
			continue;
		}
		// jump events without 2H
		if (valid == 0) {
			if (simulate_) {
				channel.taf_index = -4;
				opt.Fill();
			}
			continue;
		}
		// jump events with more than one 2H
		if (channel.taf_index > 0 && valid > 1) {
			++conflict_taf;
			if (simulate_) {
				channel.taf_index = -1;
				opt.Fill();
			}
			continue;
		}

		// reaction point
		double tx, ty, tz;
		tx = ty = tz = 0.0;
		ppac_flag = 0;
		xnum = xppac.num;
		tx = xppac.x[3];
		ty = xppac.y[3];
		tz = xppac.z[3];


		// fill extra information
		channel.entry = entry;
		t0_index[0] = t0.index[index[0]];
		t0_index[1] = t0.index[index[1]];
		csi_index = channel.taf_index < 0
			? -1
			: channel.taf_index * 2 + taf[channel.taf_index].index[0];

		// only found possible stopped 2H
		if (channel.taf_index < 0) {
			opt.Fill();
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
				t0.x[index[i]] - tx,
				t0.y[index[i]] - ty,
				t0.z[index[i]] - tz
			);
			p[i] = p[i].Unit() * momentum;
			channel.daughter_px[i] = p[i].X();
			channel.daughter_py[i] = p[i].Y();
			channel.daughter_pz[i] = p[i].Z();
			channel.daughter_r[i] = p[i].R();
			channel.daughter_theta[i] = p[i].Theta();
			channel.daughter_phi[i] = p[i].Phi();
			// other information
			channel.status[i] = t0.status[index[i]];
		}

		// rebuild real Q value
		// fill particles from TAF
		channel.recoil = 1;
		channel.recoil_charge = 1;
		channel.recoil_mass = recoil_mass_;
		channel.recoil_energy = taf[channel.taf_index].energy[0];
		channel.recoil_time = taf[channel.taf_index].time[0];
		// calculate momentum
		double recoil_momentum = MomentumFromEnergy(
			channel.recoil_energy,
			channel.recoil_charge,
			channel.recoil_mass
		);
		// recoil p
		ROOT::Math::XYZVector rp(
			taf[channel.taf_index].x[0] - tx,
			taf[channel.taf_index].y[0] - ty,
			taf[channel.taf_index].z[0] - tz
		);
		rp = rp.Unit() * recoil_momentum;
		channel.recoil_px = rp.X();
		channel.recoil_py = rp.Y();
		channel.recoil_pz = rp.Z();
		channel.recoil_r = rp.R();
		channel.recoil_theta = rp.Theta();
		channel.recoil_phi = rp.Phi();

		// fill beam from PPAC
		channel.beam_energy = 385.0;
		channel.beam_time = ppac_flag == 1 ? vppac.time[3] : xppac.time[3];
		ROOT::Math::XYZVector bp(
			ppac_flag == 1 ? vppac.px[3] : xppac.px[3],
			ppac_flag == 1 ? vppac.py[3] : xppac.py[3],
			ppac_flag == 1 ? vppac.pz[3] : xppac.pz[3]
		);
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
		q_value = channel.daughter_energy[0] + channel.daughter_energy[1]
			+ channel.recoil_energy - channel.parent_energy;

		channel.parent_px = pp.X();
		channel.parent_py = pp.Y();
		channel.parent_pz = pp.Z();
		channel.parent_r = pp.R();
		channel.parent_theta = pp.Theta();
		channel.parent_phi = pp.Phi();

		++total;
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	std::cout << "Total valid event is " << total << "\n"
		<< "Conflict TAF event " << conflict_taf << "\n"
		<< "Conflcit T0 event " << conflict_t0 << "\n";

	std::cout << "Possible 2H " << possible_2h << "\n";

	// save trees
	opt.Write();
	// fake_tree.Write();
	// close file
	t0_file.Close();
	opf.Close();
	return 0;
}



//-----------------------------------------------------------------------------
// 								C15pdChannel
//-----------------------------------------------------------------------------

C15pdChannel::C15pdChannel(unsigned int run)
: Channel(run) {
}


int C15pdChannel::Coincide() {
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
		"%s%sC15-14C-2H-%04u.root",
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
	// theta in beam coordinates
	double bct[4];
	// index of 14C in T0 event
	int t0_index;
	// coincide
	bool coincide;
	// q value
	double q_value;
	// setup output branches
	channel.SetupOutput(&opt);
	opt.Branch("beam_coord_theta", bct, "bct[num]/D");
	opt.Branch("t0_index", &t0_index, "t0_index/I");
	opt.Branch("coincide", &coincide, "co/O");
	opt.Branch("q", &q_value, "q/D");

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
		t0_index = -1;
		bool t0_valid = true;
		for (unsigned short i = 0; i < t0.num; ++i) {
			if (t0.charge[i] == 6 && t0.mass[i] == 14) {
				// find 14C
				if (t0_index < 0) {
					t0_index = i;
				} else {
					t0_valid = false;
				}
			}
		}
		if (!t0_valid && t0_index > 0) ++conflict_t0;
		// if (!t0_valid || t0_index < 0) coincide = false;
		// else coincide = true;
		if (!t0_valid || t0_index < 0) continue;
		if (t0.num == 0) coincide = false;
		else coincide = true;
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
		// if (xppac.num != 4) continue;
		// if (pow(xppac.x[3]+3.3, 2.0)+pow(xppac.y[3]-1.0, 2.0) > 225.0) {
		// 	continue;
		// }
		// if (pow(xppac.x[3]+3.0, 2.0) + pow(xppac.y[3]-1.0, 2.0) > 225.0) continue;

		// fill channel particle number
		channel.num = 2;

		if (t0_index > 0) {
			// fill particles from T0
			channel.daughter_charge[0] = t0.charge[t0_index];
			channel.daughter_mass[0] = t0.mass[t0_index];
			channel.daughter_energy[0] = t0.energy[t0_index];
			channel.daughter_time[0] = t0.time[t0_index];
			// fill p0
			ROOT::Math::XYZVector p0(
				t0.x[t0_index] - xppac.x[3],
				t0.y[t0_index] - xppac.y[3],
				t0.z[t0_index] - xppac.z[3]
			);
			double momentum = MomentumFromEnergy(
				channel.daughter_energy[0],
				channel.daughter_charge[0],
				channel.daughter_mass[0]
			);
			p0 = p0.Unit() * momentum;
			channel.daughter_px[0] = p0.X();
			channel.daughter_py[0] = p0.Y();
			channel.daughter_pz[0] = p0.Z();
			channel.daughter_r[0] = p0.R();
			channel.daughter_theta[0] = p0.Theta();
			channel.daughter_phi[0] = p0.Phi();
			// other information
			channel.status[0] = t0.status[t0_index];
		}

		// fill particles from TAF
		channel.daughter_charge[1] = taf[taf_index].charge[0];
		channel.daughter_mass[1] = taf[taf_index].mass[0];
		channel.daughter_energy[1] = taf[taf_index].energy[0];
		channel.daughter_time[1] = taf[taf_index].time[0];
		// fill p1
		ROOT::Math::XYZVector p1(
			taf[taf_index].x[0] - xppac.x[3],
			taf[taf_index].y[0] - xppac.y[3],
			taf[taf_index].z[0] - xppac.z[3]
		);
		double momentum = MomentumFromEnergy(
			channel.daughter_energy[1],
			channel.daughter_charge[1],
			channel.daughter_mass[1]
		);
		p1 = p1.Unit() * momentum;
		channel.daughter_px[1] = p1.X();
		channel.daughter_py[1] = p1.Y();
		channel.daughter_pz[1] = p1.Z();
		channel.daughter_r[1] = p1.R();
		channel.daughter_theta[1] = p1.Theta();
		channel.daughter_phi[1] = p1.Phi();

		// fill null recoil
		channel.recoil = 0;

		// fill beam from PPAC
		channel.beam_energy = 430.0;
		channel.beam_time = xppac.time[3];
		ROOT::Math::XYZVector bp(xppac.px[3], xppac.py[3], xppac.pz[3]);
		double bp_value = MomentumFromEnergy(channel.beam_energy, 6, 15);
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
		ROOT::Math::XYZVector pp = bp - p1;
		channel.parent_energy = EnergyFromMomentum(
			pp.R(), channel.parent_charge, channel.parent_mass
		);
		channel.parent_px = pp.X();
		channel.parent_py = pp.Y();
		channel.parent_pz = pp.Z();
		channel.parent_r = pp.R();
		channel.parent_theta = pp.Theta();
		channel.parent_phi = pp.Phi();

		q_value =
			430.0 + 1.0064550
			- channel.daughter_energy[1] - channel.parent_energy;

		// entry
		channel.entry = entry;
		// taf index
		channel.taf_index = taf_index;
		// fill fragments theta in beam coordinates
		bct[1] = acos(p1.Unit().Dot(bp.Unit()));


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


//-----------------------------------------------------------------------------
// 							C14ToHe4He4He6Channel
//-----------------------------------------------------------------------------

C14ToHe4He4He6Channel::C14ToHe4He4He6Channel(unsigned int run)
: Channel(run) {
}


int C14ToHe4He4He6Channel::Coincide() {
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
		"%s%sC14-4He-4He-6He-%04u.root",
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
		int valid_he6 = 0;
		unsigned short alpha_index[16];
		unsigned short he6_index[8];
		for (unsigned short i = 0; i < t0.num; ++i) {
			if (t0.charge[i] == 2){
				if (t0.mass[i] == 4) {
					alpha_index[valid_alpha] = i;
					++valid_alpha;
				} else if (t0.mass[i] == 6) {
					he6_index[valid_he6] = i;
					++valid_he6;
				}
			}
		}
		if (valid_alpha < 2 || valid_he6 < 1) continue;
		// check PPAC tracking
		if (xppac.num != 4) continue;
		if (pow(xppac.x[3]+3.3, 2.0)+pow(xppac.y[3]-1.0, 2.0) > 225.0) {
			continue;
		}

		// fill channel particle number
		channel.num = 3;
		// momentums
		std::vector<ROOT::Math::XYZVector> p;
		// fill 4He particles from T0
		for (unsigned short i = 0; i < 2; ++i) {
			channel.daughter_charge[i] = 2;
			channel.daughter_mass[i] = 4;
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
			// other information
			channel.status[i] = t0.status[alpha_index[i]];
		}

		// fill 6He particles from T0
		channel.daughter_charge[2] = 2;
		channel.daughter_mass[2] = 6;
		channel.daughter_energy[2] = t0.energy[he6_index[0]];
		channel.daughter_time[2] = t0.time[he6_index[0]];
		// calculate momentum
		double momentum = MomentumFromEnergy(
			channel.daughter_energy[2],
			channel.daughter_charge[2],
			channel.daughter_mass[2]
		);
		p.emplace_back(
			t0.x[he6_index[0]] - xppac.x[3],
			t0.y[he6_index[0]] - xppac.y[3],
			t0.z[he6_index[0]] - xppac.z[3]
		);
		p[2] = p[2].Unit() * momentum;
		channel.daughter_px[2] = p[2].X();
		channel.daughter_py[2] = p[2].Y();
		channel.daughter_pz[2] = p[2].Z();
		channel.daughter_r[2] = p[2].R();
		channel.daughter_theta[2] = p[2].Theta();
		channel.daughter_phi[2] = p[2].Phi();
		// other information
		channel.status[2] = t0.status[he6_index[0]];

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
		channel.parent_mass = 14;
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
		// parent kinetic energy
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
// 							C14ToBe10He4H1ThreeBodyChannel
//-----------------------------------------------------------------------------

C14ToBe10He4H1ThreeBodyChannel::C14ToBe10He4H1ThreeBodyChannel(unsigned int run)
: Channel(run) {
}


int C14ToBe10He4H1ThreeBodyChannel::Coincide() {
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
	unsigned short xppac_xflag;
	unsigned short xppac_yflag;
	// setup input branches
	t0.SetupInput(ipt);
	for (int i = 0; i < 6; ++i) {
		taf[i].SetupInput(ipt, "taf"+std::to_string(i)+".");
	}
	xppac.SetupInput(ipt, "xppac.");
	ipt->SetBranchAddress("xppac.xflag", &xppac_xflag);
	ipt->SetBranchAddress("xppac.yflag", &xppac_yflag);

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%sC14-10Be-4He-1H-%04u.root",
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
	int t0_index[8];
	double q_value;
	// setup output branches
	channel.SetupOutput(&opt);
	opt.Branch("t0_index", t0_index, "t0_index[num]/I");
	opt.Branch("q", &q_value, "q/D");

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
					t0.mass[i] = 10;
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
				&& (
					taf[i].mass[0] == 1
					// particle stops in ADSSD
					// taf[i].mass[0] == 0
				)
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
		// if (xppac_xflag != 0x7 || xppac_yflag != 0x7) continue;
		if (pow(xppac.x[3]+4.0, 2.0)+pow(xppac.y[3]-1.0, 2.0) > 225.0) {
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
			// other information
			channel.status[i] = t0.status[index[i]];
		}

		// rebuild real Q value
		// fill particles from TAF
		channel.recoil = 1;
		channel.recoil_charge = 1;
		channel.recoil_mass = 1;
		channel.recoil_energy = taf[taf_index].energy[0];
		channel.recoil_time = taf[taf_index].time[0];
		// calculate momentum
		double recoil_momentum = MomentumFromEnergy(
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
		rp = rp.Unit() * recoil_momentum;
		channel.recoil_px = rp.X();
		channel.recoil_py = rp.Y();
		channel.recoil_pz = rp.Z();
		channel.recoil_r = rp.R();
		channel.recoil_theta = rp.Theta();
		channel.recoil_phi = rp.Phi();

		// fill beam from PPAC
		channel.beam_energy = 385.0;
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
		q_value = channel.daughter_energy[0] + channel.daughter_energy[1]
			+ channel.recoil_energy - channel.parent_energy;

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
		// t0 index
		t0_index[0] = t0.index[index[0]];
		t0_index[1] = t0.index[index[1]];

		++total;
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	std::cout << "Total valid event is " << total << "\n"
		<< "Conflict TAF event " << conflict_taf << "\n"
		<< "Conflcit T0 event " << conflict_t0 << "\n";

	// save trees
	opt.Write();
	// close file
	t0_file.Close();
	opf.Close();
	return 0;
}


int MergeTaf(unsigned int run) {
	// TAF0 telescope file name
	TString taf0_file_name;
	taf0_file_name.Form(
		"%s%staf0-particle-ta-%04u.root",
		kGenerateDataPath,
		kParticleDir,
		run
	);
	// TAF0 telescope file
	TFile taf0_file(taf0_file_name, "read");
	// TAF0 telescope tree
	TTree *ipt = (TTree*)taf0_file.Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< taf0_file_name << " failed.\n";
		return -1;
	}
	// add TAF friends
	for (int i = 1 ; i < 6; ++i) {
		ipt->AddFriend(
			TString::Format("taf%d=tree", i),
			TString::Format(
				"%s%staf%d-particle-ta-%04u.root",
				kGenerateDataPath,
				kParticleDir,
				i,
				run
			)
		);
	}
	// input TAF particle events
	ParticleEvent taf[6];
	// setup input branches
	taf[0].SetupInput(ipt);
	for (int i = 1; i < 6; ++i) {
		taf[i].SetupInput(ipt, "taf"+std::to_string(i)+".");
	}

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%staf-particle-ta-%04u.root",
		kGenerateDataPath,
		kParticleDir,
		run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "merged TAF particles");
	// output channel event
	ParticleEvent particle;
	// setup output branches
	particle.SetupOutput(&opt);

	// total valid events
	long long total = 0;
	long long conflict_taf = 0;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Merging TAF particles   0%%");
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

		// initialize
		particle.num = 0;
		for (int i = 0; i < 6; ++i) {
			if (
				taf[i].num == 1
				&& taf[i].charge[0] == 1
				&& taf[i].mass[0] == 2
				&& taf[i].energy[0] > -9e4
			) {
				if (particle.num >= 4) continue;
				particle.charge[particle.num] = 1;
				particle.mass[particle.num] = 2;
				particle.energy[particle.num] = taf[i].energy[0];
				particle.time[particle.num] = taf[i].time[0];
				particle.x[particle.num] = taf[i].x[0];
				particle.y[particle.num] = taf[i].y[0];
				particle.z[particle.num] = taf[i].z[0];
				particle.px[particle.num] = taf[i].px[0];
				particle.py[particle.num] = taf[i].py[0];
				particle.pz[particle.num] = taf[i].pz[0];
				particle.status[particle.num] = taf[i].status[0];
				particle.index[particle.num] = i;
				++particle.num;
			}
		}
		if (particle.num == 1) {
			++total;
		} else if (particle.num > 1) {
			++conflict_taf;
		}
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	std::cout << "Total valid event is " << total << "\n"
		<< "Conflict TAF event " << conflict_taf << "\n";

	// save tree
	opt.Write();
	// close file
	taf0_file.Close();
	opf.Close();
	return 0;
}





}	// namespace ribll