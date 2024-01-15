#include <iostream>

#include <TFile.h>
#include <TGraph.h>
#include <TRandom3.h>

#include "include/defs.h"

using namespace ribll;

constexpr double kBeamEnergy = 421.0;
constexpr double kBeamEnergySigma = 3.8;
constexpr double kFragmentEnergyStep = 0.1;

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
/// @param[in] charge charge number of particle
/// @param[in] mass mass number of particle
/// @param[in] excited_energy excited energy of particle
///		default is 0.0 for ground state
/// @returns momentum value of particle
///
double MomentumFromEnergy(
	double kinetic_energy,
	unsigned short charge,
	unsigned short mass,
	double excited_energy = 0.0
) {
	// atomic mass constant
	constexpr double u = 931.494;
	double ion_mass = IonMass(charge, mass) * u + excited_energy;
	return sqrt(kinetic_energy * (kinetic_energy + 2.0*ion_mass));
}



/// @brief calculate kinetic energy from momentum considering relative effects
/// @param[in] momentum momentum of particle, in MeV
/// @param[in] mass mass number of particle
/// @returns kinetic energy of particle
///
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

int main() {
	// output file name
	TString output_file_name = TString::Format(
		"%s%stransfer-reaction.root",
		kGenerateDataPath,
		kSimulateDir
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// E-theta graph
	TGraph ground_energy_theta;
	TGraph excited_energy_theta;

	double beam_energy = kBeamEnergy;
	double beam_momentum = MomentumFromEnergy(beam_energy, 6, 15);
	double mass_energy = IonMass(1, 1) + IonMass(6, 15) - IonMass(1, 2) - IonMass(6, 14);
	double excited_energy = 0.0;
	double total_fragment_energy = beam_energy + mass_energy - excited_energy;

	for (
		double deuterium_energy = kFragmentEnergyStep;
		deuterium_energy < total_fragment_energy;
		deuterium_energy += kFragmentEnergyStep
	) {
		// get fragment energy
		double carbon_energy = total_fragment_energy - deuterium_energy;
		// calculate fragment momentum
		double deuterium_momentum = MomentumFromEnergy(deuterium_energy, 1, 2);
		double carbon_momentum = MomentumFromEnergy(carbon_energy, 6, 14);
		// cosine law
		double deuterium_cos_theta = (
			pow(deuterium_momentum, 2.0)
			+ pow(beam_momentum, 2.0)
			- pow(carbon_momentum, 2.0)
		) / (2.0 * deuterium_momentum * beam_momentum);
		if (deuterium_cos_theta > 1.0 || deuterium_cos_theta < -1.0) continue;
		double deuterium_theta = acos(deuterium_cos_theta);
		// fill to graph
		ground_energy_theta.AddPoint(deuterium_theta, deuterium_energy);
		// std::cout << deuterium_energy << " " << deuterium_theta << "  " << (
		// 		pow(deuterium_momentum, 2.0)
		// 		+ pow(carbon_momentum, 2.0)
		// 		- pow(beam_momentum, 2.0)
		// 	)
		// 	/ (2.0 * deuterium_momentum * beam_momentum) << "\n";
	}

	excited_energy = 6.093;
	total_fragment_energy = beam_energy + mass_energy - excited_energy;

	for (
		double deuterium_energy = kFragmentEnergyStep;
		deuterium_energy < total_fragment_energy;
		deuterium_energy += kFragmentEnergyStep
	) {
		// get fragment energy
		double carbon_energy = total_fragment_energy - deuterium_energy;
		// calculate fragment momentum
		double deuterium_momentum = MomentumFromEnergy(deuterium_energy, 1, 2);
		double carbon_momentum = MomentumFromEnergy(carbon_energy, 4, 14, excited_energy);
		// cosine law
		double deuterium_cos_theta = (
			pow(deuterium_momentum, 2.0)
			+ pow(beam_momentum, 2.0)
			- pow(carbon_momentum, 2.0)
		) / (2.0 * deuterium_momentum * beam_momentum);
		if (deuterium_cos_theta > 1.0 || deuterium_cos_theta < -1.0) continue;
		double deuterium_theta = acos(deuterium_cos_theta);
		// fill to graph
		excited_energy_theta.AddPoint(deuterium_theta, deuterium_energy);
	}

	// save graph
	ground_energy_theta.Write("gset");
	excited_energy_theta.Write("fset");

	// close file
	opf.Close();
	return 0;
}