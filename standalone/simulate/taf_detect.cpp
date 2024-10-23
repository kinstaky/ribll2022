#include <iostream>

#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TString.h>
#include <TGraph.h>
#include <Math/Vector3D.h>

#include "include/event/generate_event.h"
#include "include/event/dssd_event.h"
#include "include/event/csi_event.h"
#include "include/calculator/range_energy_calculator.h"

using namespace ribll;

const std::pair<double, double> tafd_phi_ranges[6] = {
	{117.6*TMath::DegToRad(), 62.4*TMath::DegToRad()},
	{57.6*TMath::DegToRad(), 2.4*TMath::DegToRad()},
	{-2.4*TMath::DegToRad(), -57.6*TMath::DegToRad()},
	{-62.4*TMath::DegToRad(), -117.6*TMath::DegToRad()},
	{-122.4*TMath::DegToRad(), -177.6*TMath::DegToRad()},
	{177.6*TMath::DegToRad(), 122.4*TMath::DegToRad()}
};
// constexpr double tafd_threshold[6] = {0.6, 0.5, 0.6, 0.6, 0.5, 0.5};
constexpr double tafd_threshold[16] = {
	0.66, 0.66, 0.7, 0.72,
	0.7, 0.72, 0.72, 0.81,
	0.84, 0.9, 0.95, 0.99,
	1.07, 1.25, 1.1, 1.1
};
// constexpr double tafd_threshold[16] = {
// 	0.6, 0.6, 0.63, 0.63,
// 	0.65, 0.65, 0.65, 0.65,
// 	0.67, 0.67, 0.67, 0.67,
// 	0.7, 0.7, 1.1, 1.1
// };
constexpr double csi_threshold[12] = {
	1200.0, 1200.0,
	1200.0, 1200.0,
	1200.0, 1400.0,
	1000.0, 1200.0,
	1100.0, 1200.0,
	1100.0, 1200.0
};

// CsI geometry parameters
// range from most distant point to z-axis (mm)
constexpr double csi_r = 174.4;
// bottom angle (alpha) of trapezoidal-shape surface
const double bottom_angle = asin(104.0 / 107.35);
// CsI range
constexpr double csi_distance = 151.8;
// CsI thickness
constexpr double csi_thickness = 31.0;

// whether CsI-9 is in low statistics like experiment
// if true, only 20% data of CsI-9 will be kept
constexpr double low_csi9 = true;


/// @brief intersection of special plane and line 特殊平面和直线的交点
/// @param[in] line_point a point in line 直线上的任意一点
/// @param[in] line_direction line direction 直线方向
/// @param[in] plane_point a point in plane 平面上任意一点
/// @param[in] plane_normal_vector normal vector of plane, z should be zero
///		平面法向量, 要求 z 方向分量为 0
/// @param[out] intersection intersection 交点
/// @returns 0 if intersection is exists, -1 otherwise
///
int IntersectionOfPlaneAndLine(
	const ROOT::Math::XYZVector &line_point,
	const ROOT::Math::XYZVector &line_direction,
	const ROOT::Math::XYZVector &plane_point,
	const ROOT::Math::XYZVector &plane_normal_vector,
	ROOT::Math::XYZVector &intersection
) {
	double k1 = plane_normal_vector.Y() / plane_normal_vector.X();
	double k2 = line_direction.X() / line_direction.Y();
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	if (fabs(plane_normal_vector.X()) < 1e-7) {
		if (
			fabs(line_direction.Y()) < 1e-7
			&& fabs(line_point.Y() - plane_point.Y()) > 1e-7
		) {
			return -1;
		}
		y = plane_point.Y();
	} else if (fabs(line_direction.Y()) < 1e-7) {
		y = line_point.Y();
	} else if (fabs(k1+k2) < 1e-7) {
		return -1;
	} else {
		y = (
			plane_point.Y() * k1
			+ line_point.Y() * k2
			+ plane_point.X()
			- line_point.Y()
		) / (k1 + k2);
	}
	if (fabs(line_direction.Y()) < 1e-7) {
		if (fabs(line_direction.X()) < 1e-7) return -1;
		if (fabs(plane_normal_vector.X() < 1e-7)) return -1;
		x = k1 * (y - plane_point.Y()) + plane_point.X();
		z = line_point.Z()
			+ (x - line_point.X()) * line_direction.Z() / line_direction.X();
	} else {
		x = line_point.X()
			+ (y - line_point.Y()) * line_direction.X() / line_direction.Y();
		z = line_point.Z()
			+ (y - line_point.Y()) * line_direction.Z() / line_direction.Y();
	}
	intersection = ROOT::Math::XYZVector(x, y, z);
	return 0;
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

	// input generate data file name
	TString generate_file_name = TString::Format(
		"%s%sgenerate-%04d.root",
		kGenerateDataPath,
		kSimulateDir,
		run
	);
	// generate data file
	TFile generate_file(generate_file_name, "read");
	// input generate tree
	TTree *generate_tree = (TTree*)generate_file.Get("tree");
	if (!generate_tree) {
		std::cerr << "Error: Get tree from "
			<< generate_file_name << " failed.\n";
		return -1;
	}
	// generate event
	GenerateEvent event;
	// setup input branches
	event.SetupInput(generate_tree);


	// output TAFD files name
	TString tafd_file_names[6];
	for (int i = 0; i < 6; ++i) {
	 	tafd_file_names[i].Form(
			"%s%stafd%d-merge-sim-ta-%04d.root",
			kGenerateDataPath,
			kMergeDir,
			i,
			run
		);
	}
	// output TAFD files
	TFile *tafd_files[6];
	for (int i = 0; i < 6; ++i) {
		tafd_files[i] = new TFile(tafd_file_names[i], "recreate");
	}
	// output TAFD trees
	TTree *tafd_trees[6];
	for (int i = 0; i < 6; ++i) {
		tafd_files[i]->cd();
		tafd_trees[i] = new TTree("tree", "simulated TAFD data");
	}
	// output TAFD merge events
	AdssdMergeEvent tafd_events[6];
	// setup output branches
	for (int i = 0; i < 6; ++i) {
		tafd_events[i].SetupOutput(tafd_trees[i]);
	}

	// output TAFCsI file name
	TString csi_file_name = TString::Format(
		"%s%stafcsi-fundamental-sim-ta-%04d.root",
		kGenerateDataPath,
		kFundamentalDir,
		run
	);
	// output TAFCsI file
	TFile csi_file(csi_file_name, "recreate");
	// output TAFCsI tree
	TTree csi_tree("tree", "simulated TAFCsI data");
	// output TAFCsI event
	CircularCsiFundamentalEvent csi_event;
	// setup output branches
	csi_event.SetupOutput(&csi_tree);

	// output TAF detect file for monitoring
	TString detect_file_name = TString::Format(
		"%s%staf-detect-%04d.root",
		kGenerateDataPath,
		kSimulateDir,
		run
	);
	// output TAF detect file
	TFile detect_file(detect_file_name, "recreate");
	// output TAF detect tree
	TTree detect_tree("tree", "TAF detect information");
	// TAF index
	int taf_index;
	// CsI index
	int csi_index;
	// particle stop status
	// 	-1: not in CsI
	//	 0: stop inside CsI
	//	 1: pass through uppper surface
	//	 2: pass through inner surface
	//	 3: pass through outer surface
	//	 4: pass through back surface
	int stop_status;
	// TAFD position
	double tafd_x, tafd_y, tafd_z;
	// TAFD strip
	int tafd_xstrip;
	// CsI front surface position
	double csi_x, csi_y, csi_z;
	// stop/pass point
	double stop_x, stop_y, stop_z;
	// range inside CsI
	double range;
	// ideal energy, should be lost in CsI if CsI is big enough
	double csi_ideal_energy;
	// CsI real energy, consider size of CsI
	double csi_real_energy;
	// setup output branches
	detect_tree.Branch("taf_index", &taf_index, "tafi/I");
	detect_tree.Branch("csi_index", &csi_index, "ci/I");
	detect_tree.Branch("tafd_x", &tafd_x, "tafdx/D");
	detect_tree.Branch("tafd_y", &tafd_y, "tafdy/D");
	detect_tree.Branch("tafd_z", &tafd_z, "tafdz/D");
	detect_tree.Branch("tafd_xstrip", &tafd_xstrip, "tafdxs/I");
	detect_tree.Branch("csi_x", &csi_x, "csix/D");
	detect_tree.Branch("csi_y", &csi_y, "csiy/D");
	detect_tree.Branch("csi_z", &csi_z, "csiz/D");
	detect_tree.Branch("stop", &stop_status, "stop/I");
	detect_tree.Branch("stop_x", &stop_x, "sx/D");
	detect_tree.Branch("stop_y", &stop_y, "sy/D");
	detect_tree.Branch("stop_z", &stop_z, "sz/D");
	detect_tree.Branch("range", &range, "range/D");
	detect_tree.Branch("csi_ideal_energy", &csi_ideal_energy, "cie/D");
	detect_tree.Branch("csi_real_energy", &csi_real_energy, "cre/D");

	// initialize random number generator
	TRandom3 generator(0);
	// initialize calculators
	elc::RangeEnergyCalculator h2_si_calculator("2H", "Si");
	elc::RangeEnergyCalculator h2_csi_calculator("2H", "CsI");

	// initialize useless variables
	for (int i = 0; i < 6; ++i) {
		for (int j = 0; j < 4; ++j) {
			tafd_events[i].time[j] = 0.0;
			tafd_events[i].decode_entry[j] = 0;
		}
	}
	csi_event.cfd_flag = 0;
	for (int i = 0; i < 12; ++i) {
		csi_event.decode_entry[i] = 0;
	}

	// pseudo random for CsI-9
	int csi_9_random = 0;

	
	// upper surface vertex
	ROOT::Math::XYZVector upper_surface_vertex[12];
	// upper surface normal vector
	ROOT::Math::XYZVector upper_surface_normal[12];
	// inner side surface (surface between two CsI) vertex
	ROOT::Math::XYZVector inner_surface_vertex[12];
	// inner side surface (surface between two CsI) normal vector
	ROOT::Math::XYZVector inner_surface_normal[12];
	// outer side surface vertex
	ROOT::Math::XYZVector outer_surface_vertex[12];
	// outer side surface normal vector
	ROOT::Math::XYZVector outer_surface_normal[12];
	// inner bottom vertex
	ROOT::Math::XYZVector inner_bottom_vertex[12];
	// CsI front trapezidal surface
	TGraph csi_front_surface[12];
	// construct plane surface for CsI
	for (size_t i = 0; i < 12; ++i) {
		// upper surface
		double upper_vertex_angle = pi * (1.0/2.0 - 1.0/3.0*(i/2));
		upper_surface_vertex[i] = ROOT::Math::XYZVector(
			csi_r * cos(upper_vertex_angle),
			csi_r * sin(upper_vertex_angle),
			csi_distance
		);
		double upper_normal_angle = i % 2 == 0
			? pi * (1.0 - (i/2)/3.0) - bottom_angle
			: bottom_angle - pi * (i/2) / 3.0;
		upper_surface_normal[i] = ROOT::Math::XYZVector(
			cos(upper_normal_angle),
			sin(upper_normal_angle),
			0.0
		);
		// inner surface
		inner_surface_vertex[i] = upper_surface_vertex[i];
		double inner_normal_angle = i % 2 == 0
			? -pi * (i/2) / 3.0
			: pi * (1.0 - (i/2)/3.0);
		inner_surface_normal[i] = ROOT::Math::XYZVector(
			cos(inner_normal_angle),
			sin(inner_normal_angle),
			0.0
		);
		// outer surface
		double outer_vertex_angle = i % 2 == 0
			? atan(64.33/157.94) + pi * (0.5 - (i/2)/3.0)
			: atan(157.94/64.33) - pi * (i/2)/3.0;
		double outer_vertex_r = sqrt(64.33*64.33 + 157.94*157.94);
		outer_surface_vertex[i] = ROOT::Math::XYZVector(
			outer_vertex_r * cos(outer_vertex_angle),
			outer_vertex_r * sin(outer_vertex_angle),
			csi_distance
		);
		double outer_normal_angle = i % 2 == 0
			? (2.0 - (i/2)/3.0) * pi - 2.0 * bottom_angle
			: 2.0 * bottom_angle - (1.0 + (i/2)/3.0) * pi;
		outer_surface_normal[i] = ROOT::Math::XYZVector(
			cos(outer_normal_angle),
			sin(outer_normal_angle),
			0.0
		);

		// bottom vertex
		inner_bottom_vertex[i] = ROOT::Math::XYZVector(
			(csi_r - 107.35) * cos(upper_vertex_angle),
			(csi_r - 107.35) * sin(upper_vertex_angle),
			csi_distance
		);
		ROOT::Math::XYZVector outer_bottom_vertex =	inner_bottom_vertex[i]
			+ (outer_surface_vertex[i] - upper_surface_vertex[i]).Unit() * 13.2;

		csi_front_surface[i].AddPoint(
			upper_surface_vertex[i].X(), upper_surface_vertex[i].Y()
		);
		csi_front_surface[i].AddPoint(
			inner_bottom_vertex[i].X(), inner_bottom_vertex[i].Y()
		);
		csi_front_surface[i].AddPoint(
			outer_bottom_vertex.X(), outer_bottom_vertex.Y()
		);
		csi_front_surface[i].AddPoint(
			outer_surface_vertex[i].X(), outer_surface_vertex[i].Y()
		);
	}

	// show start
	printf("Simulating TAF detect   0%%");
	// total entries, fow showing process
	long long entries = generate_tree->GetEntries();
	// 1/100 of entries, for showing process
	long long entry100 = entries / 100;
	// detecting simulation
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// get entry
		generate_tree->GetEntry(entry);

		// initialize
		for (int i = 0; i < 6; ++i) {
			tafd_events[i].hit = 0;
		}
		csi_event.match = false;
		for (int i = 0; i < 12; ++i) {
			csi_event.time[i] = -1e5;
			csi_event.energy[i] = 0.0;
		}

		// TAF position
		ROOT::Math::XYZVector taf_position(
			event.rx, event.ry, event.rz
		);
		tafd_x = event.rx;
		tafd_y = event.ry;
		tafd_z = event.rz;
		// TAF index
		taf_index = int(8 - event.recoil_phi / pi * 3.0) % 6;
		// middle phi
		double mid_phi = (
			tafd_phi_ranges[taf_index].second
			+ tafd_phi_ranges[taf_index].first
		) / 2.0;
		// circle center of ADSSD
		ROOT::Math::XYZVector center(
			34.4*cos(mid_phi), 34.4*sin(mid_phi), 135.0
		);
		// position relate to cirlce center
		auto relative_position = taf_position - center;
		// radius from cirlce center
		double radius = relative_position.R();
		// phi value
		double phi = relative_position.Phi();
		// out of range
		if (
			radius < 32.6
			|| radius > 135.1
			|| phi < tafd_phi_ranges[taf_index].second
			|| phi > tafd_phi_ranges[taf_index].first
		) {
			for (int i = 0; i < 6; ++i) tafd_trees[i]->Fill();
			csi_tree.Fill();
			stop_status = -1;
			detect_tree.Fill();
			continue;
		}
		// cover range of single ring strip
		double tafd_ring_width = (135.1 - 32.6) / 16.0;
		// TAFD ring strip
		int tafd_ring_strip = int((radius-32.6) / tafd_ring_width);
		tafd_xstrip = tafd_ring_strip;
		// cover range of single phi strip
		double tafd_phi_width = -55.2 * TMath::DegToRad() / 8.0;
		// TAFD phi strip
		int tafd_phi_strip = int(
			(phi - tafd_phi_ranges[taf_index].first) / tafd_phi_width
		);

		// convert strips to position back
		double merge_radius =
			tafd_ring_width * (tafd_ring_strip + 0.5) + 32.6;
		double merge_phi =
			tafd_phi_width * (tafd_phi_strip + 0.5)
			+ tafd_phi_ranges[taf_index].first;
		ROOT::Math::XYZVector merge_position(
			merge_radius * cos(merge_phi),
			merge_radius * sin(merge_phi),
			0.0
		);
		merge_position += center;

		// TAFCsI index
		csi_index = 2 * taf_index;
		if (tafd_phi_strip >= 4) csi_index += 1;

		// calculate taf energy
		double recoil_range =
			h2_si_calculator.Range(event.recoil_kinetic_after_target);
		double thickness = tafd_thickness[taf_index];
		// double thickness = 140.0;
		int taf_layer;
		double taf_lost_energy[2];
		if (recoil_range < thickness / cos(event.recoil_theta)) {
			taf_layer = 0;
			taf_lost_energy[0] = event.recoil_kinetic_after_target;
			taf_lost_energy[1] = 0.0;
		} else {
			taf_layer = 1;
			taf_lost_energy[1] = h2_si_calculator.Energy(
				recoil_range - thickness / cos(event.recoil_theta)
			);
			taf_lost_energy[0] =
				event.recoil_kinetic_after_target - taf_lost_energy[1];
		}

		// consider limited shape of CsI
		ROOT::Math::XYZVector target_point(event.target_x, event.target_y, 0.0);
		// recoil particle trace
		ROOT::Math::XYZVector recoil_direction = (taf_position - target_point).Unit();
		// calculate position on CsI front surface
		ROOT::Math::XYZVector csi_front_position = target_point
			+ recoil_direction * ((csi_distance - target_point.Z()) / recoil_direction.Z());
		csi_x = csi_front_position.X();
		csi_y = csi_front_position.Y();
		csi_z = csi_front_position.Z();
		// check recoil particle hits CsI front surface
		if (taf_layer == 0 || !csi_front_surface[csi_index].IsInside(csi_x, csi_y)) {
			stop_status = -1;
			csi_real_energy = 0.0;
			taf_layer = 0;
		} else {
			ROOT::Math::XYZVector upper_intersection, inner_intersection, outer_intersection;
			// get intersection with upper surface
			int upper_result = IntersectionOfPlaneAndLine(
				target_point, recoil_direction,
				upper_surface_vertex[csi_index], upper_surface_normal[csi_index],
				upper_intersection
			);
			int inner_result = IntersectionOfPlaneAndLine(
				target_point, recoil_direction,
				inner_surface_vertex[csi_index], inner_surface_normal[csi_index],
				inner_intersection
			);
			int outer_result = IntersectionOfPlaneAndLine(
				target_point, recoil_direction,
				outer_surface_vertex[csi_index], outer_surface_normal[csi_index],
				outer_intersection
			);

			double minz = csi_distance + csi_thickness;
			ROOT::Math::XYZVector intersection = csi_front_position
				+ recoil_direction * (csi_thickness / recoil_direction.Z());
			stop_status = 4;
			if (
				upper_result == 0
				&& upper_intersection.Z() > csi_distance
				&& upper_intersection.Z() < minz
			) {
				minz = upper_intersection.Z();
				intersection = upper_intersection;
				stop_status = 1;
			}
			if (
				inner_result == 0
				&& inner_intersection.Z() > csi_distance
				&& inner_intersection.Z() < minz
			) {
				minz = inner_intersection.Z();
				intersection = inner_intersection;
				stop_status = 2;
			}
			if (
				outer_result == 0
				&& outer_intersection.Z() > csi_distance
				&& outer_intersection.Z() < minz
			) {
				minz = outer_intersection.Z();
				intersection = outer_intersection;
				stop_status = 3;
			}

			// get range (depth) in CsI
			range = (intersection - csi_front_position).R();
			// ideal CsI energy
			csi_ideal_energy = taf_lost_energy[1];
			// ideal range
			double ideal_range = h2_csi_calculator.Range(csi_ideal_energy) / 1000.0;
			if (ideal_range > range) {
				stop_x = intersection.X();
				stop_y = intersection.Y();
				stop_z = intersection.Z();
				csi_real_energy = h2_csi_calculator.Energy(range*1000);
			} else {
				csi_real_energy = csi_ideal_energy;
				stop_status = 0;
				ROOT::Math::XYZVector stop_position =
					csi_front_position + recoil_direction * range;
				stop_x = stop_position.X();
				stop_y = stop_position.Y();
				stop_z = stop_position.Z();
			}
		}

		// consider energy resolution
		double tafd_energy = taf_lost_energy[0] + generator.Gaus(0.0, 0.05);
		double csi_energy = csi_real_energy + generator.Gaus(0.0, 1.0);
		// double tafd_energy = taf_lost_energy[0];
		// double csi_energy = taf_lost_energy[1];

		// convert energy to channel
		double csi_channel =
			power_csi_param[csi_index][0]
			* pow(csi_energy, power_csi_param[csi_index][1])
			+ power_csi_param[csi_index][2];

		// fill event
		// fill TAFD merge event
		if (tafd_energy > tafd_threshold[tafd_ring_strip]) {
			tafd_events[taf_index].hit = 1;
			tafd_events[taf_index].radius[0] = merge_position.R();
			tafd_events[taf_index].theta[0] = merge_position.Theta();
			tafd_events[taf_index].phi[0] = merge_position.Phi();
			tafd_events[taf_index].front_strip[0] = tafd_ring_strip;
			tafd_events[taf_index].back_strip[0] = tafd_phi_strip;
			tafd_events[taf_index].energy[0] = tafd_energy;
// std::cout << "Entry " << entry << ", layer " << taf_layer << ": Generate energy " << event.recoil_kinetic_after_target <<  ", range "
// 	<< recoil_range << ", taf index " << taf_index << ", csi_index " << csi_index << ", thick " << thickness
// 	<< ", effect thickness " << thickness / cos(event.recoil_theta) << ", tafd energy " << tafd_energy << ", csi energy " << csi_energy
// 	<< ", csi channel " << csi_channel << "\n";

		}
		// fill TAFCsI event
		if (taf_layer == 1 && csi_channel > csi_threshold[csi_index]) {
			if (low_csi9 && csi_index == 9) {
				if (csi_9_random == 4) {
					csi_event.match = true;
					csi_9_random = 0;
				} else {
					csi_event.match = false;
					++csi_9_random;
				}
			} else {
				csi_event.match = true;
			}
			if (csi_event.match) {
				csi_event.time[csi_index] = 0.0;
				csi_event.energy[csi_index] = csi_channel;
			}
// std::cout << "Entry " << entry << ", layer " << taf_layer << ": Generate energy " << event.recoil_kinetic_after_target <<  ", range "
// 	<< recoil_range << ", taf index " << taf_index << ", csi_index " << csi_index << ", thick " << thickness
// 	<< ", effect thickness " << thickness / cos(event.recoil_theta)  << ", tafd energy " << tafd_energy << ", csi energy " << csi_energy
// 	<< ", csi channel " << csi_channel << "\n";
		}

		// fill to tree
		for (int i = 0; i < 6; ++i) tafd_trees[i]->Fill();
		csi_tree.Fill();
		detect_tree.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");


	// save trees and close output files
	for (int i = 0; i < 6; ++i) {
		tafd_files[i]->cd();
		tafd_trees[i]->Write();
		tafd_files[i]->Close();
	}
	csi_file.cd();
	csi_tree.Write();
	csi_file.Close();
	detect_file.cd();
	for (int i = 0; i < 12; ++i) {
		csi_front_surface[i].Write(TString::Format("gf%d", i));
	}
	detect_tree.Write();
	detect_file.Close();
	// close input file
	generate_file.Close();

	return 0;
}