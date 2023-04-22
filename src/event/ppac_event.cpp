#include "include/event/ppac_event.h"

namespace ribll {

void PpacMapEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"index").c_str(), &index);
	tree->SetBranchAddress((prefix+"side").c_str(), &side);
	tree->SetBranchAddress((prefix+"cfd").c_str(), &cfd_flag);
	tree->SetBranchAddress((prefix+"time").c_str(), &time);
}


void PpacMapEvent::SetupOutput(TTree *tree) {
	tree->Branch("index", &index, "index/s");
	tree->Branch("side", &side, "side/s");
	tree->Branch("cfd", &cfd_flag, "cfd/O");
	tree->Branch("time", &time, "time/D");
}


void PpacFundamentalEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"flag").c_str(), &flag);
	tree->SetBranchAddress((prefix+"cfd").c_str(), &cfd_flag);
	tree->SetBranchAddress((prefix+"hit").c_str(), &hit);
	tree->SetBranchAddress((prefix+"xhit").c_str(), &x_hit);
	tree->SetBranchAddress((prefix+"yhit").c_str(), &y_hit);
	tree->SetBranchAddress((prefix+"x1").c_str(), x1);
	tree->SetBranchAddress((prefix+"x2").c_str(), x2);
	tree->SetBranchAddress((prefix+"y1").c_str(), y1);
	tree->SetBranchAddress((prefix+"y2").c_str(), y2);
	tree->SetBranchAddress((prefix+"anode").c_str(), anode);
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


void PpacMergeEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"xflag").c_str(), &xflag);
	tree->SetBranchAddress((prefix+"yflag").c_str(), &yflag);
	tree->SetBranchAddress((prefix+"x").c_str(), x);
	tree->SetBranchAddress((prefix+"y").c_str(), y);
}


void PpacMergeEvent::SetupOutput(TTree *tree) {
	tree->Branch("xflag", &xflag, "xflag/s");
	tree->Branch("yflag", &yflag, "yflag/s");
	tree->Branch("x", x, "x[3]/D");
	tree->Branch("y", y, "y[3]/D");
}

}