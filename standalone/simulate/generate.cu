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


__device__ void InelasticScattering(
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
}


__device__ void BreakupReaction(
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


__device__ void Rotate(
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

__device__ void VectorMultiply(double3 &vec, double multiply) {
	vec.x *= multiply;
	vec.y *= multiply;
	vec.z *= multiply;
}

__global__ void React(GenerateEvent *event, unsigned int n ) {
	unsigned int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index < n) {
		double beam_energy = beam_mass + event[index].beam_kinetic_in_target;
		double parent_energy, recoil_energy, parent_angle, recoil_angle;

		double parent_mass = beam_mass + event[index].beam_excited_energy;

		// inelastic scattering
		InelasticScattering(
			beam_mass, target_mass,
			beam_energy, event[index].beam_excited_energy,
			event[index].elastic_angle,
			parent_energy, recoil_energy, parent_angle, recoil_angle
		);
		event[index].parent_kinetic = parent_energy - parent_mass;
		event[index].recoil_kinetic_in_target = recoil_energy - target_mass;

		// rotate from beam frame to lab frame
		Rotate(
			event[index].beam_theta, event[index].beam_phi,
			parent_angle, event[index].parent_phi,
			event[index].parent_theta, event[index].parent_phi
		);
		Rotate(
			event[index].beam_theta, event[index].beam_phi,
			recoil_angle, event[index].recoil_phi,
			event[index].recoil_theta, event[index].recoil_phi
		);

		// get angle of parent relate to recoil
		double3 v_parent = {
			sin(parent_angle)*cos(event[index].parent_phi),
			sin(parent_angle)*sin(event[index].parent_phi),
			cos(parent_angle)
		};
		double v_parent_value =
			sqrt(pow(parent_energy, 2.0) - pow(parent_mass, 2.0))
			/ parent_energy;
		VectorMultiply(v_parent, v_parent_value);

		// recoil vector direction
		double3 v_recoil = {
			sin(recoil_angle)*cos(event[index].recoil_phi),
			sin(recoil_angle)*sin(event[index].recoil_phi),
			cos(recoil_angle)
		};
		double v_recoil_value =
			sqrt(pow(recoil_energy, 2.0) - pow(target_mass, 2.0))
			/ recoil_energy;
		VectorMultiply(v_recoil, v_recoil_value);

		double3 relative_v_parent_recoil = {
			v_parent.x - v_recoil.x,
			v_parent.y - v_recoil.y,
			v_parent.z - v_recoil.z
		};
		event[index].parent_recoil_angle =
			sqrt(pow(relative_v_parent_recoil.x, 2.0) + pow(relative_v_parent_recoil.y, 2.0))
			/ relative_v_parent_recoil.z;

		event[index].rz = 135.0;
		event[index].rr = event[index].rz * tan(event[index].recoil_theta);
		event[index].rx = event[index].rr * cos(event[index].recoil_phi)
			+ event[index].target_x;
		event[index].ry = event[index].rr * sin(event[index].recoil_phi)
			+ event[index].target_y;
		event[index].rr = sqrt(pow(event[index].rx, 2.0)
			+ pow(event[index].ry, 2.0));

		double fragment1_mass =
			mass_10be + event[index].fragment_excited_energy;

		double fragment1_energy, fragment2_energy;
		double fragment1_angle, fragment2_angle;
		BreakupReaction(
			parent_mass, fragment1_mass, fragment2_mass,
			parent_energy, event[index].breakup_angle,
			fragment1_energy, fragment2_energy,
			fragment1_angle, fragment2_angle
		);

		// rotate and get fragment 1 angle
		Rotate(
			event[index].parent_theta, event[index].parent_phi,
			fragment1_angle, event[index].fragment_phi_center,
			event[index].fragment_theta[0], event[index].fragment_phi[0]
		);
		// rotate and get fragment 2 angle
		Rotate(
			event[index].parent_theta, event[index].parent_phi,
			fragment2_angle, event[index].fragment_phi_center-pi,
			event[index].fragment_theta[1], event[index].fragment_phi[1]
		);
		// calculate kinetic energy for fragments
		event[index].fragment_kinetic_in_target[0] =
			fragment1_energy - fragment1_mass;
		event[index].fragment_kinetic_in_target[1] =
			fragment2_energy - fragment2_mass;


		// // get angle of parent relate to recoil
		// ROOT::Math::XYZVector v_fragment1(
		// 	sin(fragment1_angle)*cos(event[index].fragment_phi_center),
		// 	sin(fragment1_angle)*sin(event[index].fragment_phi_center),
		// 	cos(fragment1_angle)
		// );
		// v_fragment1 *=
		// 	sqrt(pow(fragment1_energy, 2.0) - pow(fragment1_mass, 2.0))
		// 	/ sqrt(pow(fragment1_energy, 2.0) + pow(fragment1_mass, 2.0));

		// ROOT::Math::XYZVector v_fragment2(
		// 	sin(fragment2_angle)*cos(event[index].fragment_phi_center-pi),
		// 	sin(fragment2_angle)*sin(event[index].fragment_phi_center-pi),
		// 	cos(fragment2_angle)
		// );
		// v_fragment2 *=
		// 	sqrt(pow(fragment2_energy, 2.0) - pow(fragment2_mass, 2.0))
		// 	/ sqrt(pow(fragment2_energy, 2.0) + pow(fragment2_mass, 2.0));

		// event[index].fragment_fragment_angle = (v_fragment1 - v_fragment2).Theta();
	}
}


__host__ void BeforeReact(
	int run,
	size_t index,
	size_t total,
	GenerateEvent &event,
	TRandom3 &generator,
	const TH1F &dist_scatter_angle,
	const TH1F &dist_breakup_angle,
	TargetEnergyCalculator &c14_target

) {
	// get reaction point x
		event.target_x = generator.Gaus(0.0, 6.0);
		// get reaction point y
		event.target_y = generator.Gaus(0.0, 6.0);
		// get beam trace at z = -800
		double beam_trace_x = generator.Gaus(0.0, 5.0);
		double beam_trace_y = generator.Gaus(0.0, 5.0);
		event.beam_phi = atan(beam_trace_y / beam_trace_x);
		if (beam_trace_x > 0) {
			event.beam_phi = beam_trace_y < 0
				? event.beam_phi + pi
				: event.beam_phi - pi;
		}
		event.beam_theta = atan(
			sqrt(pow(beam_trace_x, 2.0) + pow(beam_trace_y, 2.0)) / 800.0
		);

		// get beam kinetic energy
		event.beam_kinetic_before_target =
			generator.Gaus(389.5, 3.5);
		// reaction point depth
		event.depth = generator.Rndm();
		// consider energy loss in target
		event.beam_kinetic_in_target = c14_target.Energy(
			event.depth / cos(event.beam_theta),
			event.beam_kinetic_before_target
		);
		// get beam excited energy
		if (run == 0) {
			if (index < total / 3) {
				event.beam_excited_energy =
					12.0125 + generator.Rndm() * (20.0 - 12.0125) + 0.1;
			} else if (index < total / 3 * 2) {
				event.beam_excited_energy =
					12.0125 + 3.368
					+ generator.Rndm() * (40.0 - 12.0125 - 3.368) + 0.1;
			} else {
				event.beam_excited_energy =
					12.0125 + 6.179
					+ generator.Rndm() * (40.0 - 12.0125 - 6.179) + 0.1;
			}
		} else if (run == 1 || run == 2) {
			if (index < total / 3) {
				event.beam_excited_energy =
					12.02 + 0.2 * ((index / 10'000) % 100);
			} else if (index < total / 3 * 2) {
				event.beam_excited_energy =
					15.39 + 0.2 * ((index / 10'000) % 100);
			} else {
				event.beam_excited_energy =
					18.20 + 0.2 * ((index / 10'000) % 100);
			}
		}

		// get elastic angle theta
		event.elastic_angle = dist_scatter_angle.GetRandom();
		// get breakup angle phi
		event.parent_phi = generator.Rndm() * 2.0 * pi;
		event.recoil_phi = event.parent_phi - pi;

		// get fragment excited energy
		event.fragment_excited_energy = 0.0;
		event.fragment_state = 0;
		if (index >= total / 3 * 2) {
			event.fragment_excited_energy = 6.179;
			event.fragment_state = 2;
		} else if (index >= total / 3) {
			event.fragment_excited_energy = 3.368;
			event.fragment_state = 1;
		}
		// get breakup angle theta
		event.breakup_angle = dist_breakup_angle.GetRandom();
		// get breakup angle phi
		event.fragment_phi_center = generator.Rndm() * 2.0 * pi;
}


__host__ void AfterReact(
	GenerateEvent &event,
	TargetEnergyCalculator &he4_target,
	TargetEnergyCalculator &be10_target,
	TargetEnergyCalculator &h2_target
) {
	// consider energy loss of recoil particle in target
	event.recoil_kinetic_after_target = h2_target.Energy(
		(1.0-event.depth) / cos(event.recoil_theta),
		event.recoil_kinetic_in_target
	);

	if (event.recoil_kinetic_after_target < 0) {
		event.recoil_kinetic_after_target = 0.0;
	}

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
}

void CheckCudaError(cudaError_t error) {
    if (error != cudaSuccess) {
        std::cerr << "CUDA Error: " << cudaGetErrorString(error) << "\n";
        exit(-1);
    }
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

	constexpr size_t event_num = 3'000'000;
	// 1/100 of total event number
	size_t event100 = event_num / 100 + 1;
	GenerateEvent *generate_events = new GenerateEvent[event_num];
	// prepare events
	// show start
	printf("Preparing events   0%%");
	fflush(stdout);
	for (size_t i = 0; i < event_num; ++i) {
		// show process
		if (i % event100 == 0) {
			printf("\b\b\b\b%3lld%%", i / event100);
			fflush(stdout);
		}
		BeforeReact(
			run, i, event_num,
			generate_events[i], generator,
			dist_scatter_angle, dist_breakup_angle,
			c14_target
		);
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	GenerateEvent *cuda_events;
	// allocate memory in GPU
	CheckCudaError(cudaMalloc(&cuda_events, event_num * sizeof(GenerateEvent)));
	// copy events to GPU
	CheckCudaError(cudaMemcpy(
		cuda_events, generate_events,
		event_num * sizeof(GenerateEvent),
		cudaMemcpyHostToDevice
	));
	// simulate in GPU with CUDA
	React<<<(event_num+255)/256, 256>>>(cuda_events, event_num); 
	// copy from GPU
	CheckCudaError(cudaMemcpy(
		generate_events, cuda_events,
		event_num * sizeof(GenerateEvent),
		cudaMemcpyDeviceToHost
	));
	// free GPU memory
	CheckCudaError(cudaFree(cuda_events));

	// after reaction
	// show start
	printf("Filling events   0%%");
	fflush(stdout);
	for (size_t i = 0; i < event_num; ++i) {
		// show process
		if (i % event100 == 0) {
			printf("\b\b\b\b%3lld%%", i / event100);
			fflush(stdout);
		}
		AfterReact(
			generate_events[i],
			he4_target, be10_target, h2_target
		);
		event = generate_events[i];
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