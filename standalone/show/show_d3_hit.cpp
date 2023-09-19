#include <iostream>
#include <memory>
#include <fstream>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TCutG.h>
#include <TH1F.h>
#include <TH2F.h>

#include "include/event/t0_event.h"
#include "include/event/dssd_event.h"

using namespace ribll;

std::unique_ptr<TCutG> ReadCut(
	const char *prefix,
	const char *particle
) {
	// cut file name
	TString cut_file_name;
	cut_file_name.Form(
		"%s%scut/t0-%s-wk-b0-%s.txt",
		kGenerateDataPath,
		kParticleIdentifyDir,
		prefix,
		particle
	);
	// open cut file to read points
	std::ifstream fin(cut_file_name.Data());
	if (!fin.good()) {
		std::cerr << "Error: Open file "
			<< cut_file_name << " failed.\n";
		return nullptr;
	}
	// result
	std::unique_ptr<TCutG> result = std::make_unique<TCutG>();
	// point index
	int point;
	// point positions
	double x, y;
	// loop to read points
	while (fin.good()) {
		fin >> point >> x >> y;
		result->SetPoint(point, x, y);
	}
	// close file
	fin.close();
	return result;
}


int main() {
	// T0 file name
	TString t0_file_name = TString::Format(
		"%s%st0-telescope-ta-0618.root",
		kGenerateDataPath,
		kTelescopeDir
	);
	// T0 file
	TFile ipf(t0_file_name, "read");
	// T0 tree
	TTree *t0_tree =(TTree*)ipf.Get("tree");
	if (!t0_tree) {
		std::cerr << "Error: Get tree from "
			<< t0_file_name << " failed.\n";
		return -1;
	}
	// Add T0D3 friend
	// t0_tree->AddFriend("d3=tree", TString::Format(
	// 	"%s%st0d3-result-ta-0618.root",
	// 	kGenerateDataPath,
	// 	kNormalizeDir
	// ));
	t0_tree->AddFriend("d3=tree", TString::Format(
		"%s%st0d3-merge-ta-0618.root",
		kGenerateDataPath,
		kMergeDir
	));
	// T0 event
	T0Event t0_event;
	// T0D3 fundamental event
	// DssdFundamentalEvent d3_event;
	DssdMergeEvent d3_event;
	// setup input branches
	t0_event.SetupInput(t0_tree);
	d3_event.SetupInput(t0_tree, "d3.");

	// output file name
	TString output_file_name = TString::Format(
		"%s%sshow-d3-hit.root",
		kGenerateDataPath,
		kShowDir
	);
	TFile opf(output_file_name, "recreate");
	// histogram of T0D1D2 dE-E
	TH2F hist_d1d2_pid("d1d2s", "T0D1D2 penetrate pid", 1000, 0, 15000, 1000, 0, 15000);
	// histogram of T0D3 front-back hit
	// TH2F hist_d3_hit("d3hit", "T0D3 front back hit", 8, 0, 8, 8, 0, 8);
	// filtered t0d3 tree
	TTree t0d3_tree("tree", "4He tail in T0D1D2");
	// extra data
	double d1e;
	double d2e;
	double d1x;
	double d2x;
	double d1y;
	double d2y;
	// setup output branches
	d3_event.SetupOutput(&t0d3_tree);
	t0d3_tree.Branch("d1e", &d1e, "d1e/D");
	t0d3_tree.Branch("d2e", &d2e, "d2e/D");
	t0d3_tree.Branch("d1x", &d1x, "d1x/D");
	t0d3_tree.Branch("d2x", &d2x, "d2x/D");
	t0d3_tree.Branch("d1y", &d1y, "d1y/D");
	t0d3_tree.Branch("d2y", &d2y, "d2y/D");

	std::unique_ptr<TCutG> he4_tail = ReadCut("d1d2-tail", "He");

	for (long long entry = 0; entry < t0_tree->GetEntriesFast(); ++entry) {
		t0_tree->GetEntry(entry);
		for (unsigned short i = 0; i < t0_event.num; ++i) {
			if (t0_event.flag[i] != 0x3) continue;
			if (he4_tail->IsInside(
				t0_event.energy[i][1], t0_event.energy[i][0]
			)) {
				hist_d1d2_pid.Fill(t0_event.energy[i][1], t0_event.energy[i][0]);
				// hist_d3_hit.Fill(d3_event.back_hit, d3_event.front_hit);
				d1e = t0_event.energy[i][0];
				d2e = t0_event.energy[i][1];
				d1x = t0_event.x[i][0];
				d1y = t0_event.y[i][0];
				d2x = t0_event.x[i][1];
				d2y = t0_event.y[i][1];
				t0d3_tree.Fill();
			}
		}
	}

	// save histograms
	hist_d1d2_pid.Write();
	// hist_d3_hit.Write();
	// save tree
	t0d3_tree.Write();
	// close files
	opf.Close();
	ipf.Close();
	return 0;
}