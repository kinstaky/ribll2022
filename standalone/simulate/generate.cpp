
#include <cmath>
#include <iostream>

#include <TFile.h>
#include <TH1F.h>
#include <TRandom3.h>
#include <TString.h>
#include <TTree.h>
#include <Math/Vector3D.h>

#include "include/event/generate_event.h"
#include "include/calculator/target_energy_calculator.h"

using namespace ribll;

using elc::TargetEnergyCalculator;

constexpr double beam_mass = mass_14c;
constexpr double target_mass = mass_2h;
constexpr double fragment2_mass = mass_4he;
constexpr double target_thickness = 9.53;


/// @brief simulate process of inelastic scattering
/// @param[in] beam_mass mass of beam particles
/// @param[in] target_mass mass of target particle
/// @param[in] beam_energy full energy of beam particle
/// @param[in] excited_energy excited energy of beam particle after reaction
/// @param[in] theta elastic angle in center of mass coordinate
/// @param[out] exit_energy full energy of excited beam particle
/// @param[out] recoil_energy full energy of recoil particle
/// @param[out] exit_angle angle of excited beam particle in lab frame
/// @param[out] recoil_angle angle of recoil particle in lab frame
///
void InelasticScattering(
	double beam_mass,
	double target_mass,
	double beam_energy,
	double excited_energy,
	double theta,
	double &exit_energy,
	double &recoil_energy,
	double &exit_angle,
	double &recoil_angle
) {
	double beam_momentum = sqrt(pow(beam_energy, 2.0) - pow(beam_mass, 2.0));
	double beta_mass_center = beam_momentum / (beam_energy + target_mass);
	double gamma_mass_center =
		1.0 / sqrt(1.0 - beta_mass_center * beta_mass_center);
	double reaction_energy = sqrt(
		(beam_energy + target_mass + beam_momentum)
		* (beam_energy + target_mass - beam_momentum)
	);
	// mass of exit particle
	double exit_mass = beam_mass + excited_energy;
	// mass of recoil particle
	double recoil_mass = target_mass;
	// momentum of exit particle or recoil particle in center of mass frame
	double exit_momentum_center =
		sqrt(
			(reaction_energy - exit_mass - recoil_mass)
			* (reaction_energy - exit_mass + recoil_mass)
			* (reaction_energy + exit_mass - recoil_mass)
			* (reaction_energy + exit_mass + recoil_mass)
		) / (2.0 * reaction_energy);
	// exit momentum parallel and vertical part
	double exit_momentum_center_parallel = exit_momentum_center * cos(theta);
	double exit_momentum_center_vertical = exit_momentum_center * sin(theta);
	// exit energy in c.m.
	double exit_energy_center = sqrt(
		exit_momentum_center * exit_momentum_center
		+ exit_mass * exit_mass
	);
	// exit energy in lab frame
	exit_energy = gamma_mass_center * exit_energy_center
		+ gamma_mass_center * beta_mass_center * exit_momentum_center_parallel;
	// exit momentum parallel part in lab frame
	double exit_momentum_parallel = gamma_mass_center * exit_momentum_center_parallel
		+ gamma_mass_center * beta_mass_center * exit_energy_center;
	// exit momentum vertical part in lab frame
	double exit_momentum_vertical = exit_momentum_center_vertical;
	// calculate exit angle in lab frame
	exit_angle = fabs(atan(exit_momentum_vertical / exit_momentum_parallel));
	exit_angle = exit_momentum_parallel > 0 ?
		exit_angle : pi - exit_angle;

	// p of recoil particle in center of mass frame
	double recoil_momentum_center = exit_momentum_center;
	// recoil momentum parallel and vertical part
	double recoil_momentum_center_parallel = -recoil_momentum_center * cos(theta);
	double recoil_momentum_center_vertical = -recoil_momentum_center * sin(theta);
	// recoil energy in c.m.
	double recoil_energy_center = sqrt(
		recoil_momentum_center * recoil_momentum_center
		+ recoil_mass * recoil_mass
	);
	// recoil energy in lab frame
	recoil_energy = gamma_mass_center * recoil_energy_center
		+ gamma_mass_center * beta_mass_center * recoil_momentum_center_parallel;
	// recoil momentum parallel part in lab frame
	double recoil_momentum_parallel =
		gamma_mass_center * recoil_momentum_center_parallel
		+ gamma_mass_center * beta_mass_center * recoil_energy_center;
	// recoil momentum vertical part in lab frame
	double recoil_momentum_vertical = recoil_momentum_center_vertical;
	// recoil angle in lab frame
	recoil_angle =
		fabs(atan(recoil_momentum_vertical / recoil_momentum_parallel));
	recoil_angle = recoil_momentum_parallel > 0 ?
		recoil_angle : pi - recoil_angle;
	return;
}


/// @brief simulate breakup reaction process
/// @param[in] parent_mass mass of parent particle
/// @param[in] fragment1_mass mass of fragment particle 1
/// @param[in] fragment2_mass mass of fragment particle 2
/// @param[in] parent_energy full energy of parent particle
/// @param[in] theta breakup angle in center of mass coordinate
/// @param[out] fragment1_energy full energy of fragment particle 1
/// @param[out] fragment2_energy full energy of fragment particle 2
/// @param[out] fragment1_angle angle of fragment particle 1 in lab frame
/// @param[out] fragment2_angle angle of fragment particle 2 in lab frame
///
void BreakupReaction(
	double parent_mass,
	double fragment1_mass,
	double fragment2_mass,
	double parent_energy,
	double theta,
	double &fragment1_energy,
	double &fragment2_energy,
	double &fragment1_angle,
	double &fragment2_angle
) {
	// momentum of parent particle in lab frame
	double parent_momentum =
		sqrt(pow(parent_energy, 2.0) - pow(parent_mass, 2.0));
	// beta of center of mass
	double beta_center = parent_momentum / parent_energy;
	// gamma of center of mass
	double gamma_center = 1.0 / sqrt(1.0 - pow(beta_center, 2.0));
	// total energy in center of mass frame
	double parent_energy_center = parent_mass;
	// momentum of fragments in center of mass frame
	double fragment_momentum =
		sqrt(
			(parent_energy_center - fragment1_mass - fragment2_mass)
			* (parent_energy_center - fragment1_mass + fragment2_mass)
			* (parent_energy_center + fragment1_mass - fragment2_mass)
			* (parent_energy_center + fragment1_mass + fragment2_mass)
		) / (2.0 * parent_mass);

	// fragment1 momentum parallel part in c.m. frame
	double fragment1_momentum_center_parallel = fragment_momentum*cos(theta);
	// fragment1 momentum vertical part in c.m. frame
	double fragment1_momentum_center_vertical = fragment_momentum*sin(theta);
	// fragment1 energy in c.m. frame
	double fragment1_energy_center = sqrt(
		pow(fragment_momentum, 2.0) + pow(fragment1_mass, 2.0)
	);
	// fragment1 energy in lab frame
	fragment1_energy = gamma_center * fragment1_energy_center
		+ gamma_center * beta_center * fragment1_momentum_center_parallel;
	// fragment1 momentum parallel part in lab frame
	double fragment1_momentum_parallel =
		gamma_center * fragment1_momentum_center_parallel
		+ gamma_center * beta_center * fragment1_energy_center;
	// fragment1 momentum vertical part in lab frame
	double fragment1_momentum_vertical = fragment1_momentum_center_vertical;
	// fragment1 angle in lab frame
	fragment1_angle =
		fabs(atan(fragment1_momentum_vertical / fragment1_momentum_parallel));
	fragment1_angle = fragment1_momentum_parallel > 0 ?
		fragment1_angle : pi - fragment1_angle;

	// fragment2 momentum parallel part in c.m. frame
	double fragment2_momentum_center_parallel = -fragment_momentum*cos(theta);
	// fragment2 momentum vertical part in c.m. frame
	double fragment2_momentum_center_vertical = -fragment_momentum*sin(theta);
	// fragment2 energy in c.m. frame
	double fragment2_energy_center = sqrt(
		pow(fragment_momentum, 2.0) + pow(fragment2_mass, 2.0)
	);
	// fragment2 energy in lab frame
	fragment2_energy = gamma_center * fragment2_energy_center
		+ gamma_center * beta_center * fragment2_momentum_center_parallel;
	// fragment2 momentum parallel part in lab frame
	double fragment2_momentum_parallel =
		gamma_center * fragment2_momentum_center_parallel
		+ gamma_center * beta_center * fragment2_energy_center;
	// fragment2 momentum vertical part in lab frame
	double fragment2_momentum_vertical = fragment2_momentum_center_vertical;
	// fragment2 angle in lab frame
	fragment2_angle =
		fabs(atan(fragment2_momentum_vertical / fragment2_momentum_parallel));
	fragment2_angle = fragment2_momentum_parallel > 0 ?
		fragment2_angle : pi - fragment2_angle;

	return;
}


void Rotate(
	double parent_theta, double parent_phi,
	double fragment_theta, double fragment_phi,
	double &theta, double &phi
) {
	double x =
		sin(fragment_theta)*cos(fragment_phi)*cos(parent_theta)*cos(parent_phi)
		-sin(fragment_theta)*sin(fragment_phi)*sin(parent_phi)
		+cos(fragment_theta)*sin(parent_theta)*cos(parent_phi);
	double y =
		sin(fragment_theta)*cos(fragment_phi)*cos(parent_theta)*sin(parent_phi)
		+sin(fragment_theta)*sin(fragment_phi)*cos(parent_phi)
		+cos(fragment_theta)*sin(parent_theta)*sin(parent_phi);
	double z =
		-sin(fragment_theta)*cos(fragment_phi)*sin(parent_theta)
		+cos(fragment_theta)*cos(parent_theta);

	theta = fabs(atan(sqrt(pow(x, 2.0) + pow(y, 2.0)) / z));
	theta = z > 0 ? theta : pi - theta;
	phi = atan(y / x);
	if (x < 0) {
		phi = y > 0 ? phi + pi : phi - pi;
	}
	return;
}


int main(int argc, char **argv) {
	int run = 0;
	if (argc > 1) {
		run = atoi(argv[1]);
	}
	if (run < 0 || run > 2) {
		std::cout << "Usage: " << argv[0] << "[run]\n"
			<< "  run        run number, default is 0\n";
	}

	std::cout << "Start to simulate run " << run << "\n";

	TargetEnergyCalculator c14_target("14C", "CD2", target_thickness);
	TargetEnergyCalculator be10_target("10Be", "CD2", target_thickness);
	TargetEnergyCalculator he4_target("4He", "CD2", target_thickness);
	TargetEnergyCalculator h2_target("2H", "CD2", target_thickness);

	// generate file name
	TString generate_file_name = TString::Format(
		"%s%sgenerate-%04d.root",
		kGenerateDataPath,
		kSimulateDir,
		run
	);
	// generate file
	TFile generate_file(generate_file_name, "recreate");
	// scattering angle distribution in c.m. frame
	TH1F dist_scatter_angle("dsa", "distribution", 900, 0, pi);
	// breakup angle distribution in c.m. frame
	TH1F dist_breakup_angle("dba", "distribution", 900, 0, pi);
	// generate tree
	TTree tree("tree", "generated data");
	// generate data
	GenerateEvent event;
	// setup branches
	event.SetupOutput(&tree);

	// initialize generator
	TRandom3 generator(0);

	// parameter in scattering angle distribution
	constexpr double scatter_const = 55.0 / 180.0 * pi;
	// intialize distributions
	for (int i = 0; i < 900; ++i) {
		double x = dist_scatter_angle.GetBinCenter(i+1);
		dist_scatter_angle.SetBinContent(
			i+1, sin(x)*exp(-x/scatter_const)
		);
		dist_breakup_angle.SetBinContent(
			i+1, sin(x)
		);
	}

	double parent_energy, recoil_energy;
	double parent_angle, recoil_angle;
	double fragment1_energy, fragment2_energy;
	double fragment1_angle, fragment2_angle;


	constexpr int entries = 3'000'000;
	// 1/100 of total number of entries, for showing process
	int entry100 = entries / 100 + 1;
	// show start
	printf("Generating data   1%%");
	fflush(stdout);
	for (int entry = 0; entry < entries; ++entry) {
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3d%%", entry / entry100);
			fflush(stdout);
		}
		// get reaction point x
		event.target_x = generator.Gaus(0.0, 6.0);
		// get reaction point y
		event.target_y = generator.Gaus(0.0, 6.0);
		// get beam trace at z = -800
		double beam_trace_x = generator.Gaus(0.0, 5.0);
		double beam_trace_y = generator.Gaus(0.0, 5.0);
		event.beam_phi = atan(beam_trace_y / beam_trace_x);
		if (beam_trace_x > 0) {
			event.beam_phi = beam_trace_y < 0 ?
				event.beam_phi + pi : event.beam_phi - pi;
		}
		event.beam_theta = atan(
			sqrt(pow(beam_trace_x, 2.0) + pow(beam_trace_y, 2.0)) / 800.0
		);

		// get beam kinetic energy
		event.beam_kinetic_before_target = generator.Gaus(389.5, 3.5);
		// reaction point depth
		event.depth = generator.Rndm();
		// consider energy loss in target
		event.beam_kinetic_in_target = c14_target.Energy(
			event.depth / cos(event.beam_theta),
			event.beam_kinetic_before_target
		);
		// total energy of beam particle
		double beam_energy = beam_mass + event.beam_kinetic_in_target;
		// get beam excited energy
		if (run == 0) {
			if (entry < entries / 3) {
				event.beam_excited_energy =
					12.0125 + generator.Rndm() * (20.0 - 12.0125) + 0.1;
			} else if (entry < entries / 3 * 2) {
				event.beam_excited_energy =
					12.0125 + 3.368
					+ generator.Rndm() * (40.0 - 12.0125 - 3.368) + 0.1;
			} else {
				event.beam_excited_energy =
					12.0125 + 6.179
					+ generator.Rndm() * (40.0 - 12.0125 - 6.179) + 0.1;
			}
		} else if (run == 1 || run == 2) {
			if (entry < entries / 3) {
				event.beam_excited_energy =
					12.02 + 0.2 * ((entry / (entries / 300)) % 100);
			} else if (entry < entries / 3 * 2) {
				event.beam_excited_energy =
					15.39 + 0.2 * ((entry / (entries / 300)) % 100);
			} else {
				event.beam_excited_energy =
					18.20 + 0.2 * ((entry / (entries / 300)) % 100);
			}
		}
		double parent_mass = beam_mass + event.beam_excited_energy;

		// get elastic angle theta
		event.elastic_angle = dist_scatter_angle.GetRandom();
		// get breakup angle phi
		double parent_phi = generator.Rndm() * 2.0 * pi;
		double recoil_phi = parent_phi - pi;

		// get fragment excited energy
		event.fragment_excited_energy = 0.0;
		event.fragment_state = 0;
		if (entry >= entries / 3 * 2) {
			event.fragment_excited_energy = 6.179;
			event.fragment_state = 2;
		} else if (entry >= entries / 3) {
			event.fragment_excited_energy = 3.368;
			event.fragment_state = 1;
		}
		double fragment1_mass = mass_10be + event.fragment_excited_energy;
		// get breakup angle theta
		event.breakup_angle = dist_breakup_angle.GetRandom();
		// get breakup angle phi
		double fragment_phi_center = generator.Rndm() * 2.0 * pi;

		InelasticScattering(
			beam_mass, target_mass,
			beam_energy, event.beam_excited_energy, event.elastic_angle,
			parent_energy, recoil_energy, parent_angle, recoil_angle
		);
		event.parent_kinetic = parent_energy - parent_mass;
		event.recoil_kinetic_in_target = recoil_energy - target_mass;

		// rotate from beam frame to lab frame
		Rotate(
			event.beam_theta, event.beam_phi,
			parent_angle, parent_phi,
			event.parent_theta, event.parent_phi
		);
		Rotate(
			event.beam_theta, event.beam_phi,
			recoil_angle, recoil_phi,
			event.recoil_theta, event.recoil_phi
		);

		// get angle of parent relate to recoil
		// parent particle velocity
		ROOT::Math::XYZVector v_parent(
			sin(parent_angle)*cos(parent_phi),
			sin(parent_angle)*sin(parent_phi),
			cos(parent_angle)
		);
		v_parent *= sqrt(pow(parent_energy, 2.0) - pow(parent_mass, 2.0))
			/ parent_energy;
		// recoil particle velocity
		ROOT::Math::XYZVector v_recoil(
			sin(recoil_angle)*cos(recoil_phi),
			sin(recoil_angle)*sin(recoil_phi),
			cos(recoil_angle)
		);
		v_recoil *= sqrt(pow(recoil_energy, 2.0) - pow(target_mass, 2.0))
			/ recoil_energy;
		// angle between two velocity
		event.parent_recoil_angle = (v_parent - v_recoil).Theta();
		// beam direction
		ROOT::Math::XYZVector beam_direction(
			sin(event.beam_theta)*cos(event.beam_phi),
			sin(event.beam_theta)*sin(event.beam_phi),
			cos(event.beam_theta)
		);
		// velocity angle between relative velocity and beam
		// 母核、反冲核的相对速度矢量与入射束流方向矢量的夹角，韩家兴论文里的 θ*
		// event.angle_theta_star =
		// 	acos((v_parent - v_recoil).Unit().Dot(beam_direction));
		event.angle_theta_star = event.elastic_angle;

		// consider energy loss of recoil particle in target
		event.recoil_kinetic_after_target = h2_target.Energy(
			(1.0-event.depth) / cos(event.recoil_theta),
			event.recoil_kinetic_in_target
		);

		if (event.recoil_kinetic_after_target < 0) {
			event.recoil_kinetic_after_target = 0.0;
		}

		event.rz = 135.0;
		event.rr = event.rz * tan(event.recoil_theta);
		event.rx = event.rr * cos(event.recoil_phi) + event.target_x;
		event.ry = event.rr * sin(event.recoil_phi) + event.target_y;
		event.rr = sqrt(pow(event.rx, 2.0) + pow(event.ry, 2.0));

		BreakupReaction(
			parent_mass, fragment1_mass, fragment2_mass,
			parent_energy, event.breakup_angle,
			fragment1_energy, fragment2_energy,
			fragment1_angle, fragment2_angle
		);

		// rotate and get fragment 1 angle
		Rotate(
			event.parent_theta, event.parent_phi,
			fragment1_angle, fragment_phi_center,
			event.fragment_theta[0], event.fragment_phi[0]
		);
		// rotate and get fragment 2 angle
		Rotate(
			event.parent_theta, event.parent_phi,
			fragment2_angle, fragment_phi_center-pi,
			event.fragment_theta[1], event.fragment_phi[1]
		);
		// calculate kinetic energy for fragments
		event.fragment_kinetic_in_target[0] =
			fragment1_energy - fragment1_mass;
		event.fragment_kinetic_in_target[1] =
			fragment2_energy - fragment2_mass;


		// get angle of parent relate to recoil
		// fragment1 velocity
		ROOT::Math::XYZVector v_fragment1(
			sin(fragment1_angle)*cos(fragment_phi_center),
			sin(fragment1_angle)*sin(fragment_phi_center),
			cos(fragment1_angle)
		);
		v_fragment1 *=
			sqrt(pow(fragment1_energy, 2.0) - pow(fragment1_mass, 2.0))
			/ fragment1_energy;
		// fragment2 velocity
		ROOT::Math::XYZVector v_fragment2(
			sin(fragment2_angle)*cos(fragment_phi_center-pi),
			sin(fragment2_angle)*sin(fragment_phi_center-pi),
			cos(fragment2_angle)
		);
		v_fragment2 *=
			sqrt(pow(fragment2_energy, 2.0) - pow(fragment2_mass, 2.0))
			/ fragment2_energy;
		// angle between fragment1 and fragment2 velocity
		event.fragment_fragment_angle = (v_fragment1 - v_fragment2).Theta();
		// angle between fragments relative velocity and beam direction
		// 衰变子核的相对速度矢量和入射束流方向的夹角，韩家兴论文里的 Psi 角
		// event.angle_psi =
        //     acos((v_fragment1 - v_fragment2).Unit().Dot(beam_direction));
		event.angle_psi = event.breakup_angle;

		// 反应平面和衰变平面的夹角，即两个平面的法向量的夹角
		// event.angle_chi = acos(
		// 	(v_parent - v_recoil).Cross(beam_direction).Unit().Dot(
		// 		(v_fragment1 - v_fragment2).Cross(beam_direction).Unit()
		// 	)
		// );
		double frag_direction_theta, frag_direction_phi;
		Rotate(
			event.elastic_angle, parent_phi,
			event.breakup_angle, fragment_phi_center,
			frag_direction_theta, frag_direction_phi
		);
		event.angle_chi = fabs(frag_direction_phi - parent_phi);
		while (event.angle_chi < 0.0) event.angle_chi += pi;
		while (event.angle_chi > pi) event.angle_chi -= pi;

		// consider energy loss of fragment1 particle in target
		event.fragment_kinetic_after_target[0] = be10_target.Energy(
			(1.0-event.depth) / cos(event.fragment_theta[0]),
			event.fragment_kinetic_in_target[0]
		);
		// consider energy loss of fragment1 particle in target
		event.fragment_kinetic_after_target[1] = he4_target.Energy(
			(1.0-event.depth) / cos(event.fragment_theta[1]),
			event.fragment_kinetic_in_target[1]
		);

		// save velocity
		event.recoil_vx = v_recoil.X();
		event.recoil_vy = v_recoil.Y();
		event.recoil_vz = v_recoil.Z();
		event.parent_vx = v_parent.X();
		event.parent_vy = v_parent.Y();
		event.parent_vz = v_parent.Z();
		event.fragment_vx[0] = v_fragment1.X();
		event.fragment_vy[0] = v_fragment1.Y();
		event.fragment_vz[0] = v_fragment1.Z();
		event.fragment_vx[1] = v_fragment2.X();
		event.fragment_vy[1] = v_fragment2.Y();
		event.fragment_vz[1] = v_fragment2.Z();

		for (size_t j = 0; j < 2; ++j) {
			event.fragment_z[j] = 100.0;
			event.fragment_r[j] =
				event.fragment_z[j] * tan(event.fragment_theta[j]);
			event.fragment_x[j] =
				event.fragment_r[j] * cos(event.fragment_phi[j])
				+ event.target_x;
			event.fragment_y[j] =
				event.fragment_r[j] * sin(event.fragment_phi[j])
				+ event.target_y;
			event.fragment_r[j] = sqrt(
				pow(event.fragment_x[j], 2.0)
				+ pow(event.fragment_y[j], 2.0)
			);
		}

		// fill to tree
		tree.Fill();
	}
	// show finished
	printf("\b\b\b\b100%%\n");

	// save histograms
	dist_scatter_angle.Write();
	dist_breakup_angle.Write();
	// save tree
	tree.Write();
	// close files
	generate_file.Close();
	return 0;
}
