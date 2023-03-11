#include "include/event/ppac_event.h"

namespace ribll {

void PpacMapEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("index", &index);
	tree->SetBranchAddress("side", &side);
	tree->SetBranchAddress("cfd", &cfd_flag);
	tree->SetBranchAddress("time", &time);
}


void PpacMapEvent::SetupOutput(TTree *tree) {
	tree->Branch("index", &index, "index/s");
	tree->Branch("side", &side, "side/s");
	tree->Branch("cfd", &cfd_flag, "cfd/O");
	tree->Branch("time", &time, "time/D");
}


void PpacFundamentalEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("flag", &flag);
	tree->SetBranchAddress("cfd", &cfd_flag);
	tree->SetBranchAddress("hit", &hit);
	tree->SetBranchAddress("xhit", &x_hit);
	tree->SetBranchAddress("yhit", &y_hit);
	tree->SetBranchAddress("x1", x1);
	tree->SetBranchAddress("x2", x2);
	tree->SetBranchAddress("y1", y1);
	tree->SetBranchAddress("y2", y2);
	tree->SetBranchAddress("anode", anode);
}


void PpacFundamentalEvent::SetupOutput(TTree *tree) {
	tree->Branch("flag", &flag, "flag/i");
	tree->Branch("cfd", &cfd_flag, "cfd/s");
	tree->Branch("hit", &hit, "hit/s");
	tree->Branch("xhit", &x_hit, "xhit/s");
	tree->Branch("yhit", &y_hit, "yhit/s");
	tree->Branch("x1", x1, "x1[3]/D");
	tree->Branch("x2", x2, "x2[3]/D");
	tree->Branch("y1", y1, "y1[3]/D");
	tree->Branch("y2", y2, "y2[3]/D");
	tree->Branch("anode", anode, "a[3]/D");
}

}