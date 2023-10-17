#include "include/event/filter_event.h"

namespace ribll {

void FilterEvent::SetupInput(TTree *tree, const std::string &prefix) {
	tree->SetBranchAddress((prefix+"num").c_str(), &num);
	tree->SetBranchAddress((prefix+"pid_index").c_str(), pid_index);
	tree->SetBranchAddress((prefix+"merge_flag").c_str(), merge_flag);
	tree->SetBranchAddress(
		(prefix+"norm_front_index").c_str(), norm_front_index
	);
	tree->SetBranchAddress(
		(prefix+"norm_back_index").c_str(), norm_back_index
	);
	tree->SetBranchAddress((prefix+"front_index").c_str(), front_index);
	tree->SetBranchAddress((prefix+"back_index").c_str(), back_index);
	tree->SetBranchAddress((prefix+"front_strip").c_str(), front_strip);
	tree->SetBranchAddress((prefix+"back_strip").c_str(), back_strip);
	tree->SetBranchAddress((prefix+"front_energy").c_str(), front_energy);
	tree->SetBranchAddress((prefix+"back_energy").c_str(), back_energy);
}

void FilterEvent::SetupOutput(TTree *tree) {
	tree->Branch("num", &num, "num/I");
	tree->Branch("pid_index", pid_index, "pidi[num]/s");
	tree->Branch("merge_flag", merge_flag, "mflag[num]/s");
	tree->Branch("norm_front_index", norm_front_index, "normfi[num]/s");
	tree->Branch("norm_back_index", norm_back_index, "normbi[num]/s");
	tree->Branch("front_index", front_index, "fi[num]/s");
	tree->Branch("back_index", back_index, "bi[num]/s");
	tree->Branch("front_strip", front_strip, "fs[num]/s");
	tree->Branch("back_strip", back_strip, "bs[num]/s");
	tree->Branch("front_energy", front_energy, "fe[num]/D");
	tree->Branch("back_energy", back_energy, "be[num]/D");
}

}	// namespace ribll