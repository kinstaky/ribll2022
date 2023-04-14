#include "include/event/t0_event.h"

namespace ribll {

void T0Event::SetupInput(TTree *tree, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"num").c_str(), &num);
	tree->SetBranchAddress((prefix+"layer").c_str(), layer);
	tree->SetBranchAddress((prefix+"flag").c_str(), flag);
	tree->SetBranchAddress((prefix+"energy").c_str(), energy);
	tree->SetBranchAddress((prefix+"x").c_str(), x);
	tree->SetBranchAddress((prefix+"y").c_str(), y);
	tree->SetBranchAddress((prefix+"z").c_str(), z);
}


void T0Event::SetupOutput(TTree *tree) {
	tree->Branch("num", &num, "num/s");
	tree->Branch("layer", layer, "l[num]/s");
	tree->Branch("flag", flag, "flag[num]/s");
	tree->Branch("energy", energy, "e[num][8]/D");
	tree->Branch("x", x, "x[num][8]]/D");
	tree->Branch("y", y, "y[num][8]/D");
	tree->Branch("z", z, "z[num][8]/D");
}

}