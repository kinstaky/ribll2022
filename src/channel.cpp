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
			// other information
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
	bool hole[4];
	// input TAF particle events
	ParticleEvent taf[6];
	int sep_csi_index[6];
	// input XIA PPAC particle events
	ParticleEvent xppac;
	unsigned short xppac_xflag;
	unsigned short xppac_yflag;
	// setup input branches
	t0.SetupInput(ipt);
	ipt->SetBranchAddress("hole", hole);
	for (int i = 0; i < 6; ++i) {
		taf[i].SetupInput(ipt, "taf"+std::to_string(i)+".");
		ipt->SetBranchAddress(
			TString::Format("taf%d.csi_index", i),
			sep_csi_index + i
		);
	}
	xppac.SetupInput(ipt, "xppac.");
	ipt->SetBranchAddress("xppac.xflag", &xppac_xflag);
	ipt->SetBranchAddress("xppac.yflag", &xppac_yflag);

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
	// output fake tree
	TTree fake_tree("ftree", "channel with fake data");
	// output channel event
	ChannelEvent channel;
	int t0_index[8];
	double q_value;
	int csi_index;
	// setup output branches
	channel.SetupOutput(&opt);
	opt.Branch("t0_index", t0_index, "t0_index[num]/I");
	opt.Branch("q", &q_value, "q/D");
	opt.Branch("csi_index", &csi_index, "ci/I");
	// setup output branches for fake tree
	channel.SetupOutput(&fake_tree);
	fake_tree.Branch("t0_index", t0_index, "t0_index[num]/I");
	fake_tree.Branch("q", &q_value, "q/D");
	fake_tree.Branch("csi_index", &csi_index, "ci/I");

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
			if (hole[i]) continue;
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
					taf[i].mass[0] == 2
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
		channel.recoil_mass = 2;
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
		// csi index
		csi_index = taf_index * 2 + sep_csi_index[taf_index];
		// t0 index
		t0_index[0] = t0.index[index[0]];
		t0_index[1] = t0.index[index[1]];

		++total;
		opt.Fill();

		// rebuild fake Q value
		int rc = 12;
		int pc = 12;
		double taf_ring_width = (170.5 - 68.0) / 16.0;
		double taf_phi_width = 55.2 / 8.0;
		double tafr = sqrt(
			pow(taf[taf_index].x[0], 2.0) + pow(taf[taf_index].y[0], 2.0)
		);
		double taf_phi = atan(taf[taf_index].y[0] / taf[taf_index].x[0]);
		if (taf[taf_index].x[0] < 0) {
			taf_phi = taf[taf_index].y[0] > 0 ?
				taf_phi + 3.1415926 : taf_phi - 3.1415926;
		}
		for (int i = 0; i <= rc; ++i) {
			double fake_tafr = tafr + (i - rc/2) / double(rc) * taf_ring_width;
			for (int j = 0; j <= pc; ++j) {
				double fake_taf_phi = taf_phi
					+ (j - pc/2) / double(pc) * taf_phi_width / 180.0 * 3.1415926;
				double fake_tafx = fake_tafr * cos(fake_taf_phi);
				double fake_tafy = fake_tafr * sin(fake_taf_phi);
				// recoil p
				ROOT::Math::XYZVector fake_rp(
					fake_tafx - xppac.x[3],
					fake_tafy - xppac.y[3],
					taf[taf_index].z[0] - xppac.z[3]
				);
				fake_rp = fake_rp.Unit() * recoil_momentum;
				channel.recoil_px = fake_rp.X();
				channel.recoil_py = fake_rp.Y();
				channel.recoil_pz = fake_rp.Z();
				channel.recoil_r = fake_rp.R();
				channel.recoil_theta = fake_rp.Theta();
				channel.recoil_phi = fake_rp.Phi();

				// parent p
				ROOT::Math::XYZVector fake_pp = p[0] + p[1] + fake_rp;
				channel.parent_energy = EnergyFromMomentum(
					fake_pp.R(), channel.parent_charge, channel.parent_mass
				);
				q_value = channel.daughter_energy[0]
					+ channel.daughter_energy[1]
					+ channel.recoil_energy
					- channel.parent_energy;

				channel.parent_px = fake_pp.X();
				channel.parent_py = fake_pp.Y();
				channel.parent_pz = fake_pp.Z();
				channel.parent_r = fake_pp.R();
				channel.parent_theta = fake_pp.Theta();
				channel.parent_phi = fake_pp.Phi();

				fake_tree.Fill();
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	std::cout << "Total valid event is " << total << "\n"
		<< "Conflict TAF event " << conflict_taf << "\n"
		<< "Conflcit T0 event " << conflict_t0 << "\n";

	// save trees
	opt.Write();
	fake_tree.Write();
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
	// setup output branches
	channel.SetupOutput(&opt);
	opt.Branch("beam_coord_theta", bct, "bct[num]/D");
	opt.Branch("t0_index", &t0_index, "t0_index/I");

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
		int index = -1;
		bool t0_valid = true;
		for (unsigned short i = 0; i < t0.num; ++i) {
			if (t0.charge[i] == 6 && t0.mass[i] == 14) {
				// find 14C
				if (index < 0) {
					index = i;
				} else {
					t0_valid = false;
				}
			}
		}
		if (!t0_valid && index > 0) ++conflict_t0;
		if (!t0_valid || index < 0) continue;
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
		// if (pow(xppac.x[3]+3.3, 2.0)+pow(xppac.y[3]-1.0, 2.0) > 225.0) {
		// 	continue;
		// }
		if (pow(xppac.x[3]+3.0, 2.0) + pow(xppac.y[3]-1.0, 2.0) > 225.0) continue;

		// fill channel particle number
		channel.num = 2;
		// momentums
		std::vector<ROOT::Math::XYZVector> p;
		// fill particles from T0
		channel.daughter_charge[0] = t0.charge[index];
		channel.daughter_mass[0] = t0.mass[index];
		channel.daughter_energy[0] = t0.energy[index];
		channel.daughter_time[0] = t0.time[index];
		// fill p0
		p.emplace_back(
			t0.x[index] - xppac.x[3],
			t0.y[index] - xppac.y[3],
			t0.z[index] - xppac.z[3]
		);
		// fill particles from TAF
		channel.daughter_charge[1] = taf[taf_index].charge[0];
		channel.daughter_mass[1] = taf[taf_index].mass[0];
		channel.daughter_energy[1] = taf[taf_index].energy[0];
		channel.daughter_time[1] = taf[taf_index].time[0];
		// fill p1
		p.emplace_back(
			taf[taf_index].x[0] - xppac.x[3],
			taf[taf_index].y[0] - xppac.y[3],
			taf[taf_index].z[0] - xppac.z[3]
		);
		// calculate momentums
		for (int i = 0; i < 2; ++i) {
			double momentum = MomentumFromEnergy(
				channel.daughter_energy[i],
				channel.daughter_charge[i],
				channel.daughter_mass[i]
			);
			p[i] = p[i].Unit() * momentum;
			channel.daughter_px[i] = p[i].X();
			channel.daughter_py[i] = p[i].Y();
			channel.daughter_pz[i] = p[i].Z();
			channel.daughter_r[i] = p[i].R();
			channel.daughter_theta[i] = p[i].Theta();
			channel.daughter_phi[i] = p[i].Phi();
		}
		// other information
		channel.status[0] = t0.status[index];

		// fill null recoil
		channel.recoil = 0;

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
		channel.parent_mass = 15;
		// parent p
		ROOT::Math::XYZVector pp(0.0, 0.0, 0.0);
		for (unsigned short i = 0; i < channel.num; ++i) {
			pp += p[i];
		}
		channel.parent_energy = EnergyFromMomentum(
			pp.R(), channel.parent_charge, channel.parent_mass
		);

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
		t0_index = index;
		// fill fragments theta in beam coordinates
		for (int i = 0; i < 2; ++i) {
			bct[i] = acos(p[i].Unit().Dot(bp.Unit()));
		}


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