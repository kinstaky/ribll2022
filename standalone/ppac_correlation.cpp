#include <iostream>
#include <map>
#include <vector>

#include <TString.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

#include "include/defs.h"
#include "include/event.h"

using namespace ribll;

const long long look_window = 10000;
const long long search_window = 500;



std::vector<PPACEvent> xia_events; 

int ReadXiaEvents(int run) {
	TString file_name;
	file_name.Form("%s%sppac-corr-%04d.root", kGenerateDataPath, kCorrelationDir, run);
	TFile *ipf = new TFile(file_name, "read");
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: get tree from " << file_name << " failed.\n";
		return -1;
	}
	// input data
	PPACEvent ppac_event;
	// setup branches
	ipt->SetBranchAddress("timestamp", &ppac_event.timestamp);
	ipt->SetBranchAddress("flag", &ppac_event.flag);
	ipt->SetBranchAddress("hit", &ppac_event.hit);
	ipt->SetBranchAddress("xhit", &ppac_event.x_hit);
	ipt->SetBranchAddress("yhit", &ppac_event.y_hit);
	ipt->SetBranchAddress("x", ppac_event.x);
	ipt->SetBranchAddress("y", ppac_event.y);

	xia_events.clear();

	// show process
	printf("Reading data from xia   0%%");
	fflush(stdout);
	long long entry100 = ipt->GetEntries() / 100;
	for (long long entry = 0; entry < ipt->GetEntries(); ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		ipt->GetEntry(entry);
		xia_events.push_back(ppac_event);
	}
	printf("\b\b\b\b100%%\n");

	ipf->Close();
	return 0;
}


void PrintUsage(const char *name) {
	std::cout << "Usage: " << name << " run\n"
		<< "  run                 run number\n";
	return;
}

int main(int argc, char **argv) {
	if (argc != 2) {
		PrintUsage(argv[0]);
		return -1;
	}
	int run = atoi(argv[1]);

	ReadXiaEvents(run);

	// setup vme input file
	TString vppac_file_name;
	vppac_file_name.Form("%s%svppac-corr-%04d.root", kGenerateDataPath, kCorrelationDir, run);
	TFile *vppac_file = new TFile(vppac_file_name, "read");
	TTree *ipt = (TTree*)vppac_file->Get("tree");
	// input data
	PPACEvent vppac;
	// setup input branch
	ipt->SetBranchAddress("timestamp", &vppac.timestamp);
	ipt->SetBranchAddress("flag", &vppac.flag);
	ipt->SetBranchAddress("hit", &vppac.hit);
	ipt->SetBranchAddress("xhit", &vppac.x_hit);
	ipt->SetBranchAddress("yhit", &vppac.y_hit);
	ipt->SetBranchAddress("x", vppac.x);
	ipt->SetBranchAddress("y", vppac.y);

	// setup output file and tree
	TString align_ppac_file_name;
	align_ppac_file_name.Form("%s%sppac-align-%04d.root", kGenerateDataPath, kAlignDir, run);
	TFile *opf = new TFile(align_ppac_file_name, "recreate");
	TH1F *hist_look_window = new TH1F("ht", "look window", 1000, -look_window, look_window);
	TTree *opt = new TTree("tree", "tree of align ppac");
	// output data
	PPACEvent xppac;
	// setup output branches
	opt->Branch("xia_timestamp", &xppac.timestamp, "xts/L");
	opt->Branch("vme_timestamp", &vppac.timestamp, "vts/L");
	opt->Branch("xia_flag", &xppac.flag, "xflag/I");
	opt->Branch("vme_flag", &vppac.flag, "vflag/I");
	opt->Branch("xia_hit", &xppac.hit, "xhit/s");
	opt->Branch("vme_hit", &vppac.hit, "vhit/s");
	opt->Branch("xia_xhit", &xppac.x_hit, "xxhit/s");
	opt->Branch("vme_xhit", &vppac.x_hit, "vxhit/s");
	opt->Branch("xia_yhit", &xppac.y_hit, "xyhit/s");
	opt->Branch("vme_yhit", &vppac.y_hit, "vyhit/s");
	opt->Branch("xia_x", xppac.x, TString::Format("xx[%llu]/D", ppac_num).Data());
	opt->Branch("vme_x", vppac.x, TString::Format("vx[%llu]/D", ppac_num).Data());
	opt->Branch("xia_y", xppac.y, TString::Format("xy[%llu]/D", ppac_num).Data());
	opt->Branch("vme_y", vppac.y, TString::Format("vy[%llu]/D", ppac_num).Data());


	// show process
	printf("Correlating ppac from xia and vme   0%%");
	fflush(stdout);
	long long entry100 = ipt->GetEntries() / 100;
	// correlation rate
	long long correlation_num = 0;
	for (long long entry = 0; entry < ipt->GetEntries(); ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}
		ipt->GetEntry(entry);

		int local_hit = 0;
		for (
			auto iter = std::lower_bound(
				xia_events.begin(),
				xia_events.end(),
				vppac.timestamp - look_window,
				[](const PPACEvent &x, const long long &ts) {
					return x.timestamp < ts;
				}
			);
			iter != xia_events.end();
			++iter
		) {
			// stop
			if (iter->timestamp >= vppac.timestamp + look_window) break;
			// record look window
			hist_look_window->Fill(vppac.timestamp - iter->timestamp);
			// jump if out of search
			if (iter->timestamp <= vppac.timestamp - search_window) continue;
			if (iter->timestamp >= vppac.timestamp + search_window) continue;
			
			++local_hit;
			xppac = *iter;
		}

		if (local_hit == 1) {
			opt->Fill();
			++correlation_num;
		}
	}
	printf("\b\b\b\b100%%\n");
	hist_look_window->Write();
	opt->Write();
	opf->Close();

	std::cout << "vme ppac correlation rate " << correlation_num << " / " << ipt->GetEntries()
		<< " " << double(correlation_num) / double(ipt->GetEntries()) << "\n";
	
	vppac_file->Close();
	return 0;
}