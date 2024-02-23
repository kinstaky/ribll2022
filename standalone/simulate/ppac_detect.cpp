#include <iostream>

#include <TFile.h>
#include <TString.h>
#include <TRandom3.h>

#include "include/event/generate_event.h"
#include "include/event/ppac_event.h"
#include "include/event/tof_event.h"

using namespace ribll;

constexpr double xppac_x_efficiency = 0.8;
constexpr double xppac_y_efficiency = 0.85;
constexpr double vppac_efficiency = 0.95;
constexpr double vme_correct_x[3] = {-0.25, 1.5, 0.25};
constexpr double vme_correct_y[3] = {-0.75, 0.0, -0.75};
constexpr double xia_offset[2][3] = {{-0.20, 0.11, 9.80}, {9.43, -14.4, 1.79}};
constexpr double vme_sum_mean[2][3] = {{57.0, 64.0, 68.0}, {57.5, 66.0, 68.0}};
constexpr double vme_sum_sigma = 1.05;
// first index for x(0) and y(1),
// second index for x or y index,
// third index for anode index
constexpr double xia_sum_mean[2][3][3] = {
	// {x0-a0, x0-a1, x0-a2}, {x1-a0, x1-a1, x1-a2}, {x2-a0, x2-a1, x2-a2}
	{{107.0, 108.0, 80.0}, {103.0, 103.0, 75.0}, {114.0, 115.0, 86.0}},
	// {y0-a0, y0-a1, y0-a2}, {y1-a0, y1-a1, y1-a2}, {y2-a0, y2-a1, y2-a2}
	{{106.0, 107.0, 78.0}, {110.0, 111.0, 82.0}, {122.0, 122.0, 94.0}}
};
constexpr double xia_sum_sigma[2][3][3] = {
	// {x0-a0, x0-a1, x0-a2}, {x1-a0, x1-a1, x1-a2}, {x2-a0, x2-a1, x2-a2}
	{{0.58, 0.85, 0.86}, {0.85, 0.64, 0.88}, {0.85, 0.86, 0.60}},
	// {y0-a0, y0-a1, y0-a2}, {y1-a0, y1-a1, y1-a2}, {y2-a0, y2-a1, y2-a2}
	{{0.55, 0.83, 0.85}, {0.83, 0.59, 0.86}, {0.82, 0.83, 0.57}}
};
constexpr double xia_anode_mean[3] = {-482.5, -482.8, -468.6};
constexpr double xia_anode_sigma[3] = {0.43, 0.44, 0.44};


int main(int argc, char **argv) {
	bool vppac = false;
	if (argc == 1) {
		vppac = false;
	} else if (argc == 2 && argv[1][0] == '-' && argv[1][1] == 'v') {
		vppac = true;
	} else {
		std::cout << "Usage: " << argv[0] << " [options]\n"
			<< "Options:\n"
			<< "  -v        VME PPAC detect, default is XIA PPAC.\n";
		return -1;
	}

	// input generate data file name
	TString generate_file_name = TString::Format(
		"%s%sgenerate.root",
		kGenerateDataPath,
		kSimulateDir
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

	// VME TOF evenet
	TofFundamentalEvent vtof_event;
	if (vppac) {
		// VME TOF file name
		TString vtof_file_name = TString::Format(
			"%s%svtof-fundamental-sim-ta-0000.root",
			kGenerateDataPath,
			kFundamentalDir
		);
		// add friend
		generate_tree->AddFriend("vtof=tree", vtof_file_name);
		// setup input event
		vtof_event.SetupInput(generate_tree, "vtof.");
	}

	// output file name
	TString output_file_name = TString::Format(
		"%s%s%cppac-fundamental-sim-ta-0000.root",
		kGenerateDataPath,
		kFundamentalDir,
		vppac ? 'v' : 'x'
	);
	// output file
	TFile opf(output_file_name, "recreate");
	// output tree
	TTree opt("tree", "simulated PPAC data");
	// output event
	PpacFundamentalEvent ppac;
	// merge x and y for debug
	double merge_x[3], merge_y[3];
	// setup output branches
	ppac.SetupOutput(&opt);
	// opt.Branch("merge_x", merge_x, "mx[3]/D");
	// opt.Branch("merge_y", merge_y, "my[3]/D");

	// initialize random number generator
	TRandom3 generator(0);

	// show start
	printf("Simulating PPAC detect   0%%");
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

		// initialize output event
		ppac.flag = ppac.cfd_flag = 0;
		ppac.hit = ppac.x_hit = ppac.y_hit = 0;

		// PPAC position
		for (int i = 0; i < 3; ++i) {
			// X position
			double ppacx =
				ppac_xz[i] * tan(event.beam_theta) * cos(event.beam_phi)
				+ event.target_x;
			// consdier strips
			ppacx = int(ppacx + 50.5) - 50.0 + generator.Gaus(0.0, 0.35);
			// Y position
			double ppacy =
				ppac_yz[i] * tan(event.beam_theta) * sin(event.beam_phi)
				+ event.target_y;
			// consdier strips
			ppacy = int(ppacy + 50.5) - 50.0 + generator.Gaus(0.0, 0.35);

			// efficiency check
			for (int j = 0; j < 5; ++j) {
				double random = generator.Rndm();
				// failed to pass efficiency check
				if (
					(vppac && random > vppac_efficiency)
					|| (!vppac && j < 2 && random > xppac_x_efficiency)
					|| (!vppac && j >= 2 && j < 4 && random > xppac_y_efficiency)
					|| (!vppac && j == 4 && random > 0.95)
				) continue;
				// pass efficiency check
				ppac.flag |= 0x1 << (i*5+j);
				++ppac.hit;
			}

			// back to merge level
			merge_x[i] = ppacx * 4.0 + vme_correct_x[i];
			merge_y[i] = ppacy * -4.0 + vme_correct_y[i];
			// XIA system offset
			if (!vppac) {
				merge_x[i] -= xia_offset[0][i];
				merge_y[i] -= xia_offset[1][i];
			}

			if (vppac) {
				// get sum
				double sumx = generator.Gaus(vme_sum_mean[0][i], vme_sum_sigma);
				double sumy = generator.Gaus(vme_sum_mean[1][i], vme_sum_sigma);

				// X direction fundamental event
				double tmp_x = sumx + vtof_event.time[1] + vtof_event.time[1];
				ppac.x1[i] = (tmp_x + merge_x[i]) / 2.0;
				ppac.x2[i] = (tmp_x - merge_x[i]) / 2.0;

				// Y direction fundamental event
				double tmp_y = sumy + vtof_event.time[1] + vtof_event.time[1];
				ppac.y1[i] = (tmp_y + merge_y[i]) / 2.0;
				ppac.y2[i] = (tmp_y - merge_y[i]) / 2.0;
			} else {
				// get sum
				double sumx = generator.Gaus(
					xia_sum_mean[0][i][i], xia_sum_sigma[0][i][i]
				);
				double sumy = generator.Gaus(
					xia_sum_mean[1][i][i], xia_sum_sigma[1][i][i]
				);

				// get anode time
				double anode = generator.Gaus(xia_anode_mean[i], xia_anode_sigma[i]);

				// X directon fundamentaal event
				double tmp_x = sumx + anode + anode;
				ppac.x1[i] = (tmp_x + merge_x[i]) / 2.0;
				ppac.x2[i] = (tmp_x - merge_x[i]) / 2.0;

				// Y direction fundamental event
				double tmp_y = sumy + anode + anode;
				ppac.y1[i] = (tmp_y + merge_y[i]) / 2.0;
				ppac.y2[i] = (tmp_y - merge_y[i]) / 2.0;

				// anode
				ppac.anode[i] = anode;
			}
		}

		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save tree
	opt.Write();
	// close file
	opf.Close();
	generate_file.Close();

	return 0;
}