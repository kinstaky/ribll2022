#include "include/event/dssd_event.h"

namespace ribll {

void DssdMapEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"side").c_str(), &side);
	tree->SetBranchAddress((prefix+"strip").c_str(), &strip);
	tree->SetBranchAddress((prefix+"cfd").c_str(), &cfd_flag);
	tree->SetBranchAddress((prefix+"time").c_str(), &time);
	tree->SetBranchAddress((prefix+"energy").c_str(), &energy);
	tree->SetBranchAddress((prefix+"decode_entry").c_str(), &decode_entry);
}


void DssdMapEvent::SetupOutput(TTree *tree) {
	tree->Branch("side", &side, "side/s");
	tree->Branch("strip", &strip, "s/s");
	tree->Branch("cfd", &cfd_flag, "cfd/O");
	tree->Branch("time", &time, "t/D");
	tree->Branch("energy", &energy, "e/D");
	tree->Branch("decode_entry", &decode_entry, "decode_entry/L");
}


void DssdFundamentalEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"front_hit").c_str(), &front_hit);
	tree->SetBranchAddress((prefix+"back_hit").c_str(), &back_hit);
	tree->SetBranchAddress((prefix+"cfd").c_str(), &cfd_flag);
	tree->SetBranchAddress((prefix+"front_strip").c_str(), front_strip);
	tree->SetBranchAddress((prefix+"back_strip").c_str(), back_strip);
	tree->SetBranchAddress((prefix+"front_time").c_str(), front_time);
	tree->SetBranchAddress((prefix+"back_time").c_str(), back_time);
	tree->SetBranchAddress((prefix+"front_energy").c_str(), front_energy);
	tree->SetBranchAddress((prefix+"back_energy").c_str(), back_energy);
	// tree->SetBranchAddress(
	// 	(prefix+"front_decode_entry").c_str(), front_decode_entry
	// );
	// tree->SetBranchAddress(
	// 	(prefix+"back_decode_entry").c_str(), back_decode_entry
	// );
}


void DssdFundamentalEvent::SetupOutput(TTree *tree) {
	tree->Branch("front_hit", &front_hit, "fhit/s");
	tree->Branch("back_hit", &back_hit, "bhit/s");
	tree->Branch("cfd", &cfd_flag, "cfd/s");
	tree->Branch("front_strip", front_strip, "fs[fhit]/s");
	tree->Branch("back_strip", back_strip, "bs[bhit]/s");
	tree->Branch("front_time", front_time, "ft[fhit]/D");
	tree->Branch("back_time", back_time, "bt[bhit]/D");
	tree->Branch("front_energy", front_energy, "fe[fhit]/D");
	tree->Branch("back_energy", back_energy, "be[bhit]/D");
	// tree->Branch("front_decode_entry", front_decode_entry, "fde[fhit]/L");
	// tree->Branch("back_decode_entry", back_decode_entry, "bde[bhit]/L");
}


void DssdTimeEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"front_hit").c_str(), &front_hit);
	tree->SetBranchAddress((prefix+"back_hit").c_str(), &back_hit);
	tree->SetBranchAddress(
		(prefix+"front_time_flag").c_str(), front_time_flag
	);
	tree->SetBranchAddress(
		(prefix+"back_time_flag").c_str(), back_time_flag
	);
}


void DssdTimeEvent::SetupOutput(TTree *tree) {
	tree->Branch("front_hit", &front_hit, "fhit/s");
	tree->Branch("back_hit", &back_hit, "bhit/s");
	tree->Branch("front_time_flag", front_time_flag, "ftf[fhit]/I");
	tree->Branch("back_time_flag", back_time_flag, "btf[bhit]/I");
}


void DssdNormalizeEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"front_hit").c_str(), &front_hit);
	tree->SetBranchAddress((prefix+"back_hit").c_str(), &back_hit);
	tree->SetBranchAddress((prefix+"front_strip").c_str(), front_strip);
	tree->SetBranchAddress((prefix+"back_strip").c_str(), back_strip);
	tree->SetBranchAddress((prefix+"front_energy").c_str(), front_energy);
	tree->SetBranchAddress((prefix+"back_energy").c_str(), back_energy);
}


void DssdNormalizeEvent::SetupOutput(TTree *tree) {
	tree->Branch("front_hit", &front_hit, "fhit/s");
	tree->Branch("back_hit", &back_hit, "bhit/s");
	tree->Branch("front_strip", front_strip, "fs[fhit]/s");
	tree->Branch("back_strip", back_strip, "bs[bhit]/s");
	tree->Branch("front_energy", front_energy, "fe[fhit]/D");
	tree->Branch("back_energy", back_energy, "be[bhit]/D");
}


void DssdMergeEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"hit").c_str(), &hit);
	tree->SetBranchAddress((prefix+"case").c_str(), &case_tag);
	tree->SetBranchAddress((prefix+"x").c_str(), x);
	tree->SetBranchAddress((prefix+"y").c_str(), y);
	tree->SetBranchAddress((prefix+"z").c_str(), z);
	tree->SetBranchAddress((prefix+"energy").c_str(), energy);
	tree->SetBranchAddress((prefix+"time_flag").c_str(), time_flag);
}


void DssdMergeEvent::SetupOutput(TTree *tree) {
	tree->Branch("hit", &hit, "hit/s");
	tree->Branch("case", &case_tag, "case/i");
	tree->Branch("x", x, "x[hit]/D");
	tree->Branch("y", y, "y[hit]/D");
	tree->Branch("z", z, "z[hit]/D");
	tree->Branch("energy", energy, "e[hit]/D");
	tree->Branch("time_flag", time_flag, "tf[hit]/I");
}


void AdssdMergeEvent::SetupInput(
	TTree *tree,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"hit").c_str(), &hit);
	tree->SetBranchAddress((prefix+"radius").c_str(), radius);
	tree->SetBranchAddress((prefix+"theta").c_str(), theta);
	tree->SetBranchAddress((prefix+"phi").c_str(), phi);
	tree->SetBranchAddress((prefix+"energy").c_str(), energy);
	tree->SetBranchAddress((prefix+"decode_entry").c_str(), decode_entry);
}


void AdssdMergeEvent::SetupOutput(TTree *tree) {
	tree->Branch("hit", &hit, "hit/s");
	tree->Branch("radius", radius, "r[hit]/D");
	tree->Branch("theta", theta, "theta[hit]/D");
	tree->Branch("phi", phi, "phi[hit]/D");
	tree->Branch("energy", energy, "e[hit]/D");
	tree->Branch("decode_entry", decode_entry, "de[hit]/L");
}


}	// namespace ribll