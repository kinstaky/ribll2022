#include <iostream>
#include <iomanip>

#include <Math/Vector3D.h>
#include <TDirectoryFile.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>

#include "include/event/t0_event.h"
#include "include/event/particle_type_event.h"

using namespace ribll;

const ROOT::Math::XYZVector d1_center{0.0, 0.0, 100.0};
const ROOT::Math::XYZVector d2_center{0.0, 0.0, 111.76};
const ROOT::Math::XYZVector d3_center{0.0, 0.0, 123.52};

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " run\n"
			<< "  run               Set the run number.\n";
		return -1;
	}
	unsigned int run = atoi(argv[1]);

	// t0 file name
	TString t0_file_name;
	t0_file_name.Form(
		"%s%st0-telescope-ta-%04u.root",
		kGenerateDataPath, kTelescopeDir, run
	);
	// d1 file
	TFile t0_file(t0_file_name, "read");
	// d1 tree
	TTree *tree = (TTree*)t0_file.Get("tree");
	if (!tree) {
		std::cerr << "Error: Get tree from "
			<< t0_file_name << " failed.\n";
		return -1;
	}
	// particle type file name
	TString type_file_name;
	type_file_name.Form(
		"%s%st0-particle-type-ta-%04u.root",
		kGenerateDataPath, kParticleIdentifyDir, run
	);
	tree->AddFriend("type=tree", type_file_name);
	// input t0 event
	T0Event t0_event;
	// input particle type event
	ParticleTypeEvent type_event;
	// setup branches
	t0_event.SetupInput(tree);
	type_event.SetupInput(tree, "type.");

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%st0-center-%04u.root",
		kGenerateDataPath, kShowDir, run
	);
	// output file
	TFile opf(output_file_name, "recreate");

	// lower bound to search
	double lower_bound = -10.0;
	// step of distance
	double step = 1.0;
	// only one extreme value in cos(theta) VS distance
	bool one_extreme_value = true;
	// loop layer
	int loop = 0;
	while (one_extreme_value && loop < 10) {
		// directory name
		TString dir_name;
		dir_name.Form("d%i", loop);
		// create new subdirectory in root file
		TDirectoryFile dir(dir_name, dir_name);
		dir.cd();
		// histograms of cos(theta) in different d2 distance
		std::vector<TH1F> hist_cos_theta;
		// average cos(theta) VS distance
		TGraph cos_theta_vs_distance;

		// maximum value of cos(theta)
		double max_cos_theta = -1.0;
		// distance when get the maximum cos(theta)
		double max_cos_distance = 0.0;
		// last cos(theta)
		double last_cos_theta = -1.0;
		// whether have reached the maximum value of cos(theta)
		bool reach_max_cos = false;
		// loop the d2 distance
		for (size_t i = 0; i < 20; ++i) {
			double distance = lower_bound + i*step;
			hist_cos_theta.emplace_back(
				TString::Format("h%ld", i), "cos(#theta)",
				100, 0.999, 1
			);
			// average cos(theta)
			double average_cos_theta = 0.0;
			// total valid events
			long long count = 0;
			// loop events and fill the histogram
			for (long long entry = 0; entry < tree->GetEntriesFast(); ++entry) {
				tree->GetEntry(entry);
				if (
					t0_event.num < 1 || t0_event.flag[0] != 0x7
					|| type_event.mass[0] <= 0 || type_event.charge[0] <= 0
				) continue;
				// d1 particle position
				ROOT::Math::XYZVector d1_pos(
					t0_event.x[0][0], t0_event.y[0][0], d1_center.Z()
				);
				// d2 particle position
				ROOT::Math::XYZVector d2_pos(
					t0_event.x[0][1]-1.45738+distance, t0_event.y[0][1]-0.31811, d2_center.Z()
				);
				// d3 particle position
				// ROOT::Math::XYZVector d3_pos(
				// 	t0_event.x[0][2]-2.02242, t0_event.y[0][2]-0.03429, d3_center.Z()
				// );
				double cos_theta = d1_pos.Dot(d2_pos) / (d1_pos.R() * d2_pos.R());
				// double cos_theta = d1_pos.Dot(d3_pos) / (d1_pos.R() * d3_pos.R());
				// double cos_theta = d2_pos.Dot(d3_pos) / (d2_pos.R() * d3_pos.R());
				hist_cos_theta[i].Fill(cos_theta);
				average_cos_theta += cos_theta;
				++count;
			}
			average_cos_theta /= double(count);
			cos_theta_vs_distance.AddPoint(distance, average_cos_theta);
			if (average_cos_theta > max_cos_theta) {
				max_cos_theta = average_cos_theta;
				max_cos_distance = distance;
			}
			if (average_cos_theta < last_cos_theta && !reach_max_cos) {
				reach_max_cos = true;
			} else if (average_cos_theta > last_cos_theta && reach_max_cos) {
				one_extreme_value = false;
			}
			last_cos_theta = average_cos_theta;
			// print average cos_theta
			// std::cout << "d2 distance " << d2_distance << " mm, cos(theta) "
			// 	<< std::setprecision(20) << average_cos_theta << "\n";
		}
		if (one_extreme_value) {
			std::cout << "Max cos(theta) at " << std::setprecision(12)
				<< max_cos_distance << "\n";
		}
		// save histograms
		dir.cd();
		for (auto &hist : hist_cos_theta) {
			hist.Write();
		}
		// save graph
		cos_theta_vs_distance.Write("gcos");
		dir.Write();

		// update for next iteration
		lower_bound = max_cos_distance - step;
		step *= 0.1;
		++loop;
		hist_cos_theta.clear();
		cos_theta_vs_distance.Set(0);
	}
	// close files
	opf.Close();
	t0_file.Close();
	return 0;
}
