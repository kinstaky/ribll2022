#include "include/event/ppac_event.h"

namespace ribll {

void PpacMapEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("index", &index);
	tree->SetBranchAddress("side", &side);
	tree->SetBranchAddress("time", &time);
}


void PpacMapEvent::SetupOutput(TTree *tree) {
	tree->Branch("index", &index, "index/s");
	tree->Branch("side", &side, "side/s");
	tree->Branch("time", &time, "time/D");
}


void PpacFundamentalEvent::SetupInput(TTree *tree) {
	tree->SetBranchAddress("flag", &flag);
	tree->SetBranchAddress("hit", &hit);
	tree->SetBranchAddress("xhit", &x_hit);
	tree->SetBranchAddress("yhit", &y_hit);
	tree->SetBranchAddress("x", x);
	tree->SetBranchAddress("y", y);
}


void PpacFundamentalEvent::SetupOutput(TTree *tree) {
	tree->Branch("flag", &flag, "flag/i");
	tree->Branch("hit", &hit, "hit/s");
	tree->Branch("xhit", &x_hit, "xhit/s");
	tree->Branch("yhit", &y_hit, "yhit/s");
	tree->Branch("x", x, TString::Format("x[%llu]/D", ppac_num));
	tree->Branch("y", y, TString::Format("y[%llu]/D", ppac_num));
}

}