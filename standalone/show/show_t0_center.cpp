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

#include "include/event/dssd_event.h"

using namespace ribll;

const ROOT::Math::XYZVector d1_center{0.0, 0.0, 10.0};

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " run\n"
			<< "  run               Set the run number.\n";
		return -1;
	}
	unsigned int run = atoi(argv[1]);

	// t0d1 file name
	TString d1_file_name;
	d1_file_name.Form(
		"%s%st0d1-merge-ta-%04u.root",
		kGenerateDataPath, kMergeDir, run
	);
	// d1 file
	TFile d1_file(d1_file_name, "read");
	// d1 tree
	TTree *tree = (TTree*)d1_file.Get("tree");
	if (!tree) {
		std::cerr << "Error: Get tree from "
			<< d1_file_name << " failed.\n";
		return -1;
	}
	// t0d2 file name
	TString d2_file_name;
	d2_file_name.Form(
		"%s%st0d2-merge-ta-%04u.root",
		kGenerateDataPath, kMergeDir, run
	);
	tree->AddFriend("d2=tree", d2_file_name);
	// input d1 event
	DssdMergeEvent d1_event;
	// input d2 event
	DssdMergeEvent d2_event;
	// setup branches
	d1_event.SetupInput(tree);
	d2_event.SetupInput(tree, "d2.");

	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%st0-center-%04u.root",
		kGenerateDataPath, kShowDir, run
	);
	// output file
	TFile opf(output_file_name, "recreate");

	// lower bound to search
	double lower_bound = 110.0;
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
			double d2_distance = lower_bound + i*step;
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
				if (d1_event.hit != 1 || d2_event.hit != 1) continue;
				// d1 particle position
				ROOT::Math::XYZVector d1_pos(
					d1_event.x[0]+1.0, d1_event.y[0]+1.0, d1_event.z[0]
				);
				// d2 particle position
				ROOT::Math::XYZVector d2_pos(
					d2_event.x[0], d2_event.y[0], d2_distance
				);
				double cos_theta = d1_pos.Dot(d2_pos) / (d1_pos.R() * d2_pos.R());
				hist_cos_theta[i].Fill(cos_theta);
				average_cos_theta += cos_theta;
				++count;
			}
			average_cos_theta /= double(count);
			cos_theta_vs_distance.AddPoint(d2_distance, average_cos_theta);
			if (average_cos_theta > max_cos_theta) {
				max_cos_theta = average_cos_theta;
				max_cos_distance = d2_distance;
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
	d1_file.Close();
	return 0;
}
