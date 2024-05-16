#include <iostream>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TMarker.h>
#include <TPaveText.h>
#include <TString.h>
#include <TTree.h>

#include "include/event/dssd_event.h"

using namespace ribll;

int main(int argc, char **argv) {
	if (argc != 3) {
		std::cout << "Usage: " << argv[0] << " run end_run\n"
			<< "  run            start of run\n"
			<< "  end_run         end of run, inclusive\n";
		return -1;
	}
	int run = atoi(argv[1]);
	int end_run = atoi(argv[2]);

	// output file name
	TString output_file_name = TString::Format(
		"%s%sshow-d2-shift-%04d-%04d.root",
		kGenerateDataPath,
		kShowDir,
		run,
		end_run
	);
	// output file
	TFile opf(output_file_name, "recreate");
	TGraph graph_front[32];
	std::vector<TGraph> graph_front_run[32];
	for (int fs = 0; fs < 32; ++fs) {
		graph_front_run[fs].resize(end_run-run+1);
	}

	for (int i = run; i <= end_run; ++i) {
		if (i == 628) continue;
		// input file name
		TString input_file_name = TString::Format(
			"%s%st0d2-fundamental-ta-%04d.root",
			kGenerateDataPath,
			kFundamentalDir,
			i
		);
		// input file
		TFile ipf(input_file_name, "read");
		// input tree
		TTree *ipt = (TTree*)ipf.Get("tree");
		if (!ipt) {
			std::cerr << "Error: Get tree from "
				<< input_file_name << " failed.\n";
			continue;
		}
		// input data
		DssdFundamentalEvent event;
		// setup input branches
		event.SetupInput(ipt);

		// show process
		std::cout << "Processing run " << i << "\n";
		for (long long entry = 0; entry < ipt->GetEntriesFast(); ++entry) {
			ipt->GetEntry(entry);
			if (event.front_hit != 1 || event.back_hit != 1) continue;
			if (event.back_strip[0] != 19) continue;
			if (event.front_strip[0] < 10 || event.front_strip[0] >= 20) continue;
			if (event.front_energy[0] > 20000) continue;
			if (event.back_energy[0] > 20000) continue;
			graph_front[event.front_strip[0]].AddPoint(
				event.front_energy[0], event.back_energy[0]
			);
			for (int j = run; j <= end_run; ++j) {
				if (j == 628) continue;
				graph_front_run[event.front_strip[0]][j-run].AddPoint(
					event.front_energy[0], event.back_energy[0]
				);
			}
			TMarker *marker = new TMarker(
				event.front_energy[0], event.back_energy[0], 20
			);
			marker->SetMarkerColor(kRed);
			graph_front_run[event.front_strip[0]][i-run]
				.GetListOfFunctions()->Add(marker);
		}
		// close input file
		ipf.Close();
	}

	// save gif
	for (int i = 10; i < 20; ++i) {
		std::cout << "Storing front strip " << i << "\n";
		TCanvas c1;
		for (int j = run; j <= end_run; ++j) {
			if (j == 628) continue;
			graph_front_run[i][j-run].Draw("AP");
			TPaveText *pt = new TPaveText(17'000, 2'000, 20'000, 4'000);
			pt->AddText(TString::Format("%d", j));
			pt->SetTextColor(kRed);
			pt->SetBorderSize(0);
			pt->SetFillColor(0);
			pt->SetFillStyle(0);
			pt->Draw();
			if (j == run) {
				c1.Print(TString::Format(
					"%s%st0d2-shift-fs%d.gif",
					kGenerateDataPath,
					kShowDir,
					i
				));
			} else {
				c1.Print(TString::Format(
					"%s%st0d2-shift-fs%d.gif+50",
					kGenerateDataPath,
					kShowDir,
					i
				));
			}
		}
		c1.Print(TString::Format(
			"%s%st0d2-shift-fs%d.gif++",
			kGenerateDataPath,
			kShowDir,
			i
		));
	}


	opf.cd();
	// save graphs
	for (int i = 10; i < 20; ++i) {
		graph_front[i].Write(TString::Format("gf%d", i));
	}
	for (int i = 10; i < 20; ++i) {
		for (int j = run; j <= end_run; ++j) {
			if (j == 628) continue;
			graph_front_run[i][j-run].Write(TString::Format(
				"gf%d_%d", i, j
			));
		}
	}
	// close file
	opf.Close();
	return 0;
}