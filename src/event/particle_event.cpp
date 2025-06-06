#include "include/event/particle_event.h"

namespace ribll {

void ParticleEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"num").c_str(), &num);
	tree->SetBranchAddress((prefix+"charge").c_str(), charge);
	tree->SetBranchAddress((prefix+"mass").c_str(), mass);
	tree->SetBranchAddress((prefix+"layer").c_str(), layer);
	tree->SetBranchAddress((prefix+"energy").c_str(), energy);
	tree->SetBranchAddress((prefix+"time").c_str(), time);
	tree->SetBranchAddress((prefix+"x").c_str(), x);
	tree->SetBranchAddress((prefix+"y").c_str(), y);
	tree->SetBranchAddress((prefix+"z").c_str(), z);
	tree->SetBranchAddress((prefix+"px").c_str(), px);
	tree->SetBranchAddress((prefix+"py").c_str(), py);
	tree->SetBranchAddress((prefix+"pz").c_str(), pz);
	tree->SetBranchAddress((prefix+"status").c_str(), status);
	tree->SetBranchAddress((prefix+"index").c_str(), index);
	tree->SetBranchAddress((prefix+"straight").c_str(), straight);
}


void ParticleEvent::SetupOutput(TTree *tree) {
	tree->Branch("num", &num, "num/I");
	tree->Branch("charge", charge, "Z[num]/s");
	tree->Branch("mass", mass, "A[num]/s");
	tree->Branch("layer", layer, "layer[num]/s");
	tree->Branch("energy", energy, "e[num]/D");
	tree->Branch("time", time, "t[num]/D");
	tree->Branch("x", x, "x[num]/D");
	tree->Branch("y", y, "y[num]/D");
	tree->Branch("z", z, "z[num]/D");
	tree->Branch("px", px, "px[num]/D");
	tree->Branch("py", py, "py[num]/D");
	tree->Branch("pz", pz, "pz[num]/D");
	tree->Branch("status", status, "status[num]/I");
	tree->Branch("index", index, "index[num]/I");
	tree->Branch("straight", straight, "straight[num]/O");
}

}	// namespace ribll