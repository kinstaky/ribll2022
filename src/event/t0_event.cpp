#include "include/event/t0_event.h"

namespace ribll {

void T0Event::SetupInput(TTree *tree, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"num").c_str(), &num);
	tree->SetBranchAddress((prefix+"layer").c_str(), layer);
	tree->SetBranchAddress((prefix+"flag").c_str(), flag);
	tree->SetBranchAddress((prefix+"energy").c_str(), energy);
	tree->SetBranchAddress((prefix+"time").c_str(), time);
	tree->SetBranchAddress((prefix+"x").c_str(), x);
	tree->SetBranchAddress((prefix+"y").c_str(), y);
	tree->SetBranchAddress((prefix+"z").c_str(), z);
	tree->SetBranchAddress((prefix+"sflag").c_str(), &ssd_flag);
	tree->SetBranchAddress((prefix+"ssd_energy").c_str(), ssd_energy);
	tree->SetBranchAddress((prefix+"status").c_str(), status);
	tree->SetBranchAddress((prefix+"points").c_str(), points);
}


void T0Event::SetupOutput(TTree *tree) {
	tree->Branch("num", &num, "num/s");
	tree->Branch("layer", layer, "l[num]/s");
	tree->Branch("flag", flag, "flag[num]/s");
	tree->Branch("energy", energy, "e[num][4]/D");
	tree->Branch("time", time, "t[num][4]/D");
	tree->Branch("x", x, "x[num][4]]/D");
	tree->Branch("y", y, "y[num][4]/D");
	tree->Branch("z", z, "z[num][4]/D");
	tree->Branch("sflag", &ssd_flag, "sflag/s");
	tree->Branch("ssd_energy", ssd_energy, "se[3]/D");
	tree->Branch("status", status, "status[num]/I");
	tree->Branch("points", points, "points[num]/I");
}

}