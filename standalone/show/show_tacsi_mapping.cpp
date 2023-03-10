#include <iostream>
#include <string>

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TH1F.h>
#include <TPaveText.h>
#include <TCanvas.h>

#include "include/defs.h"

using namespace ribll;

// main canvas
TCanvas *c1;
// output file name
TString output_file_name;

int Show(unsigned int run, const std::string &name) {
	// input file name
	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-fundamental-%04u.root",
		kGenerateDataPath, kFundamentalDir, name.c_str(), run
	);
	TFile *ipf = new TFile(input_file_name, "read");
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// csi file name
	TString csi_file_name;
	csi_file_name.Form(
		"%s%s%s-fundamental-%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		(name.substr(0, 3) + "csi").c_str(),
		run
	);
	if (!ipt->AddFriend("ct=tree", csi_file_name)) {
		std::cerr << "Error: Add friend from "
			<< csi_file_name << " failed.\n";
		return -1;
	}

	// adssd front hit, input data
	unsigned short fhit;
	// adssd back hit, input data
	unsigned short bhit;
	// adssd back strip, input data
	unsigned short bs[8];
	// csi time, input data
	double csi_time[12];
	// setup input branches
	ipt->SetBranchAddress("front_hit", &fhit);
	ipt->SetBranchAddress("back_hit", &bhit);
	ipt->SetBranchAddress("ct.time", csi_time);
	ipt->SetBranchAddress("back_strip", bs);

	// histograms to store coincident entries with different CsI(Tl)
	TH1F *h[2];
	h[0] = new TH1F("hls", (name + " low strips").c_str(), 12, 0, 12);
	h[1] = new TH1F("hhs", (name + " high strips").c_str(), 12, 0, 12);

	// total entries
	long long entries = ipt->GetEntries();
	// loop to get coincident entries
	for (long long entry = 0; entry < entries; ++entry) {
		ipt->GetEntry(entry);
		if (fhit != 1 || bhit != 1) continue;
		for (int i = 0; i < 12; ++i) {
			if (csi_time[i] > 0) {
				h[bs[0]>=4]->Fill(i);
			}
		}
	}

	// print to output pdf file
	h[0]->Draw();
	c1->Print(output_file_name, ("Title:" + name + " low strips").c_str());
	h[1]->Draw();
	c1->Print(output_file_name, ("Title:" + name + " high strips").c_str());

	return 0;
}


int main(int argc, char **argv) {
	// check arguments
	if (argc != 2) {
		std::cout << "Usage: " << argv[0] << " run\n"
			<< "  run                  Run number.\n";
		return -1;
	}
	unsigned int run = atoi(argv[1]);

	// initialize global variables
	c1 = new TCanvas;
	output_file_name.Form(
		"%s%s/tacsi-mapping-%04u.pdf",
		kGenerateDataPath, kShowDir, run
	);

	// add first title page to output pdf file
	TPaveText *first_page = new TPaveText(0.05, 0.1, 0.95, 0.9);
	first_page->AddText("PDF for showing TAF and TAB csi mapping");
	first_page->Draw();
	c1->Print(output_file_name + "(", "Title:start");

	// show TAF
	for (int i = 0; i < 6; ++i) {
		std::string detector = "taf" + std::to_string(i);
		if (Show(run, detector)) {
			std::cerr << "Error: Show " << detector << " failed.\n";
		}
	}

	// show TAB
	for (int i = 0; i < 6; ++i) {
		std::string detector = "tab" + std::to_string(i);
		if (Show(run, detector)) {
			std::cerr << "Error: Show " << detector << " failed.\n";
		}
	}

	// add last ending page
	TCanvas *c2 = new TCanvas;
	TPaveText *last_page = new TPaveText(0.05, 0.1, 0.95, 0.9);
	last_page->AddText("END");
	last_page->Draw();
	c2->Print(output_file_name + ")", "Title:end");

	return 0;
}