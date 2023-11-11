#include "include/event/t0_event.h"

namespace ribll {

void T0Event::SetupInput(TTree *tree, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"num").c_str(), &num);
	tree->SetBranchAddress((prefix+"layer").c_str(), layer);
	tree->SetBranchAddress((prefix+"flag").c_str(), flag);
	tree->SetBranchAddress((prefix+"charge").c_str(), charge);
	tree->SetBranchAddress((prefix+"mass").c_str(), mass);
	tree->SetBranchAddress((prefix+"energy").c_str(), energy);
	tree->SetBranchAddress((prefix+"time").c_str(), time);
	tree->SetBranchAddress((prefix+"x").c_str(), x);
	tree->SetBranchAddress((prefix+"y").c_str(), y);
	tree->SetBranchAddress((prefix+"z").c_str(), z);
	tree->SetBranchAddress((prefix+"sflag").c_str(), &ssd_flag);
	tree->SetBranchAddress((prefix+"ssd_energy").c_str(), ssd_energy);
	tree->SetBranchAddress((prefix+"status").c_str(), status);
	tree->SetBranchAddress((prefix+"points").c_str(), points);
	tree->SetBranchAddress((prefix+"dssd_flag").c_str(), dssd_flag);
	tree->SetBranchAddress((prefix+"hole").c_str(), hole);
}


void T0Event::SetupOutput(TTree *tree) {
	tree->Branch("num", &num, "num/I");
	tree->Branch("layer", layer, "l[num]/S");
	tree->Branch("flag", flag, "flag[num]/s");
	tree->Branch("charge", charge, "Z[num]/s");
	tree->Branch("mass", mass, "A[num]/s");
	tree->Branch("energy", energy, "e[num][3]/D");
	tree->Branch("time", time, "t[num][3]/D");
	tree->Branch("x", x, "x[num][3]/D");
	tree->Branch("y", y, "y[num][3]/D");
	tree->Branch("z", z, "z[num][3]/D");
	tree->Branch("sflag", &ssd_flag, "sflag/s");
	tree->Branch("ssd_energy", ssd_energy, "se[3]/D");
	tree->Branch("status", status, "status[num]/I");
	tree->Branch("points", points, "points[num]/I");
	tree->Branch("dssd_flag", dssd_flag, "dssd_flag[num][3]/s");
	tree->Branch("hole", hole, "hole[num]/O");
}

}