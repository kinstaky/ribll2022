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
	// accurate mass in MeV/C^2
	double ion_mass = IonMass(charge, mass) * u;

	return sqrt(momentum*momentum + ion_mass*ion_mass) - ion_mass;
}


//-----------------------------------------------------------------------------
// 								Channel
//-----------------------------------------------------------------------------

Channel::Channel(
	unsigned int run,
	const std::vector<std::string> &particles
)
: run_(run)
, particles_(particles){

	for (const auto &p : particles) {
		if (p.size() < 2 || p.size() > 4) {
			throw std::runtime_error("Invalid nuclear " + p);
		}
		if (p[0] < '0' || p[0] > '9') {
			throw std::runtime_error("Invalid nuclear" + p);
		}
		if (p.size() == 2) {
			masses_.push_back(p[0]-'0');
			auto search = element.find(p.substr(1, 1));
			if (search == element.end()) {
				throw std::runtime_error("Invalid nuclear" + p);
			}
			charges_.push_back(search->second);
		} else if (p.size() == 3) {
			if (p[1] < '0' || p[1] > '9') {
				masses_.push_back(p[0] - '0');
				auto search = element.find(p.substr(1, 2));
				if (search == element.end()) {
					throw std::runtime_error("Invalid nuclear" + p);
				}
				charges_.push_back(search->second);
			} else {
				unsigned short mass = (p[0] - '0') * 10 + (p[1] - '0');
				masses_.push_back(mass);
				auto search = element.find(p.substr(2, 1));
				if (search == element.end()) {
					throw std::runtime_error("Invalid nuclear" + p);
				}
				charges_.push_back(search->second);
			}
		} else if (p.size() == 4) {
			if (p[1] < '0' || p[1] > '9') {
				throw std::runtime_error("Invalid nuclear" + p);
			}
			unsigned short mass = (p[0] - '0') * 10 + (p[1] - '0');
			masses_.push_back(mass);
			auto search = element.find(p.substr(2, 2));
			if (search == element.end()) {
				throw std::runtime_error("Invalid nuclear" + p);
			}
			charges_.push_back(search->second);
		}
	}
}

int Channel::Coincide() {
	std::cerr << "Error: Channel::Coincide not implemented yet.\n";
	return -1;
}

//-----------------------------------------------------------------------------
// 								T0TAFChannel
//-----------------------------------------------------------------------------

T0TAFChannel::T0TAFChannel(
	unsigned int run,
	const std::vector<std::string> &particles
)
: Channel(run, particles) {

}

constexpr bool recoil = false;
int T0TAFChannel::Coincide() {
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
	// input T0 particle events
	ParticleEvent t0;
	// input TAF particle events
	ParticleEvent taf[6];
	// setup input branches
	t0.SetupInput(ipt);
	for (int i = 0; i < 6; ++i) {
		taf[i].SetupInput(ipt, "taf"+std::to_string(i)+".");
	}

	// output file name
	std::string particle_names;
	for (const auto &p : particles_) particle_names += p + "-";
	TString output_file_name;
	output_file_name.Form(
		"%s%st0taf-%s%04u.root",
		kGenerateDataPath,
		kChannelDir,
		particle_names.c_str(),
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

	// particles number in T0
	unsigned short t0_num = particles_.size() - 1;

	// total valid events
	long long total = 0;
	long long conflict_taf = 0;

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

		if (t0.num != t0_num) continue;
		bool t0_valid = true;
		for (unsigned short i = 0; i < t0_num; ++i) {
			if (t0.charge[i] != charges_[i] || t0.mass[i] != masses_[i]) {
				t0_valid = false;
				break;
			}
		}
		if (!t0_valid) continue;
		// check TAF
		int taf_index = -1;
		int valid = 0;
		for (int i = 0; i < 6; ++i) {
			if (
				taf[i].num == 1
				&& taf[i].charge[0] == charges_[t0_num]
				&& taf[i].mass[0] == masses_[t0_num]
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

		// momentums
		std::vector<ROOT::Math::XYZVector> p;
		// fill particles from T0
		for (unsigned short i = 0; i < t0_num; ++i) {
			channel.charge[i] = t0.charge[i];
			channel.mass[i] = t0.mass[i];
			channel.energy[i] = t0.energy[i];
			// calculate momentum
			double momentum = MomentumFromEnergy(
				channel.energy[i],
				channel.charge[i],
				channel.mass[i]
			);
			p.emplace_back(t0.x[i], t0.y[i], t0.z[i]);
			p[i] = p[i].Unit() * momentum;
		}
		// fill particles from TAF
		channel.charge[t0_num] = taf[taf_index].charge[0];
		channel.mass[t0_num] = taf[taf_index].mass[0];
		channel.energy[t0_num] = taf[taf_index].energy[0];
		// calculate momentum
		double momentum = MomentumFromEnergy(
			channel.energy[t0_num],
			channel.charge[t0_num],
			channel.mass[t0_num]
		);
		p.emplace_back(
			taf[taf_index].x[0], taf[taf_index].y[0], taf[taf_index].z[0]
		);
		p[t0_num] = p[t0_num].Unit() * momentum;
		// fill parent nuclear
		channel.charge[particles_.size()] = 0;
		channel.mass[particles_.size()] = 0;
		if (recoil) {
			for (size_t i = 0; i < particles_.size()-1; ++i) {
				channel.charge[particles_.size()] += channel.charge[i];
				channel.mass[particles_.size()] += channel.mass[i];
			}
		} else {
			for (size_t i = 0; i < particles_.size(); ++i) {
				channel.charge[particles_.size()] += channel.charge[i];
				channel.mass[particles_.size()] += channel.mass[i];
			}
		}
		p.emplace_back(0.0, 0.0, 0.0);
		for (size_t i = 0; i < particles_.size(); ++i) {
			p[particles_.size()] += p[i];
		}
		channel.energy[particles_.size()] =
			EnergyFromMomentum(
				p[particles_.size()].R(),
				channel.charge[particles_.size()],
				channel.mass[particles_.size()]
			);

		// fill momentum
		for (size_t i = 0; i <= particles_.size(); ++i) {
			channel.px[i] = p[i].X();
			channel.py[i] = p[i].Y();
			channel.pz[i] = p[i].Z();
			channel.r[i] = p[i].R();
			channel.theta[i] = p[i].Theta();
			channel.phi[i] = p[i].Phi();
		}
		// fill channel particle number
		channel.num = particles_.size() + 1;

		++total;
		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	std::cout << "Total valid event is " << total << "\n"
		<< "Conflict TAF event " << conflict_taf << "\n";

	// save tree
	opt.Write();
	// close file
	t0_file.Close();
	opf.Close();
	return 0;
}


//-----------------------------------------------------------------------------
// 								T0Channel
//-----------------------------------------------------------------------------

T0Channel::T0Channel(
	unsigned int run,
	const std::vector<std::string> &particles
)
: Channel(run, particles) {
}

int T0Channel::Coincide() {
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
	// input T0 particle events
	ParticleEvent t0;
	//setup input branches
	t0.SetupInput(ipt);

	// output file name
	std::string particle_names;
	for (const auto &p : particles_) particle_names += p + "-";
	TString output_file_name;
	output_file_name.Form(
		"%s%st0-%s%04u.root",
		kGenerateDataPath,
		kChannelDir,
		particle_names.c_str(),
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
		if (t0.num != particles_.size()) continue;
		bool t0_valid = true;
		for (unsigned short i = 0; i < particles_.size(); ++i) {
			if (t0.charge[i] != charges_[i] || t0.mass[i] != masses_[i]) {
				t0_valid = false;
				break;
			}
		}
		if (!t0_valid) continue;

		// momentums
		std::vector<ROOT::Math::XYZVector> p;
		// fill particles from T0
		for (unsigned short i = 0; i < particles_.size(); ++i) {
			channel.charge[i] = t0.charge[i];
			channel.mass[i] = t0.mass[i];
			channel.energy[i] = t0.energy[i];
			// calculate momentum
			double momentum = MomentumFromEnergy(
				channel.energy[i],
				channel.charge[i],
				channel.mass[i]
			);
			p.emplace_back(t0.x[i], t0.y[i], t0.z[i]);
			p[i] = p[i].Unit() * momentum;
		}
		// fill parent nuclear
		channel.charge[particles_.size()] =
			channel.charge[0] + channel.charge[1];
		channel.mass[particles_.size()] = channel.mass[0] + channel.mass[1];
		p.emplace_back(0.0, 0.0, 0.0);
		for (size_t i = 0; i < particles_.size(); ++i) {
			p[particles_.size()] += p[i];
		}
		channel.energy[particles_.size()] =
			EnergyFromMomentum(
				p[particles_.size()].R(),
				channel.charge[particles_.size()],
				channel.mass[particles_.size()]
			);

		// fill momentum
		for (size_t i = 0; i <= particles_.size(); ++i) {
			channel.px[i] = p[i].X();
			channel.py[i] = p[i].Y();
			channel.pz[i] = p[i].Z();
			channel.r[i] = p[i].R();
			channel.theta[i] = p[i].Theta();
			channel.phi[i] = p[i].Phi();
		}
		// fill particle numebr
		channel.num = particles_.size() + 1;

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

}	// namespace ribll