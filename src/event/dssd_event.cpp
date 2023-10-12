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
	tree->SetBranchAddress(
		(prefix+"front_decode_entry").c_str(), front_decode_entry
	);
	tree->SetBranchAddress(
		(prefix+"back_decode_entry").c_str(), back_decode_entry
	);
	tree->SetBranchAddress(
		(prefix+"front_fundamental_index").c_str(), front_fundamental_index
	);
	tree->SetBranchAddress(
		(prefix+"back_fundamental_index").c_str(), back_fundamental_index
	);
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
	tree->Branch("front_decode_entry", front_decode_entry, "fde[fhit]/L");
	tree->Branch("back_decode_entry", back_decode_entry, "bde[bhit]/L");
	tree->Branch(
		"front_fundamental_index", front_fundamental_index, "ffi[fhit]/s"
	);
	tree->Branch(
		"back_fundamental_index", back_fundamental_index, "bfi[bhit]/s"
	);
}


void DssdFundamentalEvent::Sort() {
	SortSide(0, 0);
	SortSide(1, 0);
}


void DssdFundamentalEvent::SortSide(size_t side, unsigned short index) {
	if (side == 0) {
		// front side
		if (index == front_hit) return;
		double max_energy = front_energy[index];
		unsigned short max_index = index;
		for (unsigned short i = index+1; i < front_hit; ++i) {
			if (front_energy[i] > max_energy) {
				max_energy = front_energy[i];
				max_index = i;
			}
		}
		Swap(side, index, max_index);
		SortSide(side, index+1);
	} else {
		// back side
		if (index == back_hit) return;
		double max_energy = back_energy[index];
		unsigned short max_index = index;
		for (unsigned short i = index+1; i < back_hit; ++i) {
			if (back_energy[i] > max_energy) {
				max_energy = back_energy[i];
				max_index = i;
			}
		}
		Swap(side, index, max_index);
		SortSide(side, index+1);
	}
}


void DssdFundamentalEvent::Swap(
	size_t side,
	unsigned short i,
	unsigned short j
) {
	if (side == 0) {
		// front side
		// swap cfd flag
		unsigned short t1 = (cfd_flag & (1 << i)) >> i;
		unsigned short t2 = (cfd_flag & (1 << j)) >> j;
		cfd_flag &= ~(1 << i);
		cfd_flag &= ~(1 << j);
		cfd_flag |= t1 << j;
		cfd_flag |= t2 << i;
		// swap strip
		unsigned short ts = front_strip[i];
		front_strip[i] = front_strip[j];
		front_strip[j] = ts;
		// swap time
		double tt = front_time[i];
		front_time[i] = front_time[j];
		front_time[j] = tt;
		// swap energy
		double te = front_energy[i];
		front_energy[i] = front_energy[j];
		front_energy[j] = te;
		// swap decode entry
		long long tde = front_decode_entry[i];
		front_decode_entry[i] = front_decode_entry[j];
		front_decode_entry[j] = tde;
		// swap fundamental index
		unsigned short tfi = front_fundamental_index[i];
		front_fundamental_index[i] = front_fundamental_index[j];
		front_fundamental_index[j] = tfi;
	} else {
		// back side
		// swap cfd flag
		unsigned short t1 = (cfd_flag & (1 << (i+8))) >> (i+8);
		unsigned short t2 = (cfd_flag & (1 << (j+8))) >> (j+8);
		cfd_flag &= ~(1 << (i+8));
		cfd_flag &= ~(1 << (j+8));
		cfd_flag |= t1 << (j+8);
		cfd_flag |= t2 << (i+8);
		// swap strip
		unsigned short ts = back_strip[i];
		back_strip[i] = back_strip[j];
		back_strip[j] = ts;
		// swap time
		double tt = back_time[i];
		back_time[i] = back_time[j];
		back_time[j] = tt;
		// swap energy
		double te = back_energy[i];
		back_energy[i] = back_energy[j];
		back_energy[j] = te;
		// swap decode entry
		long long tde = back_decode_entry[i];
		back_decode_entry[i] = back_decode_entry[j];
		back_decode_entry[j] = tde;
		// swap fundamental index
		unsigned short tbi = back_fundamental_index[i];
		back_fundamental_index[i] = back_fundamental_index[j];
		back_fundamental_index[j] = tbi;
	}
}


void DssdFundamentalEvent::Erase(size_t side, size_t index) {
	// move front events
	if (side == 0) {
		for (unsigned short i = index+1; i < front_hit; ++i) {
			front_strip[i-1] = front_strip[i];
			front_energy[i-1] = front_energy[i];
			front_time[i-1] = front_time[i];
			cfd_flag |= (cfd_flag >> 1) & (1<<(i-1));
		}
		--front_hit;
	} else {
		for (unsigned short i = index+1; i < back_hit; ++i) {
			back_strip[i-1] = back_strip[i];
			back_energy[i-1] = back_energy[i];
			back_time[i-1] = back_time[i];
			cfd_flag |= (cfd_flag >> 1) & (1<<(i+8-1));
		}
		--back_hit;
	}
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
	tree->SetBranchAddress((prefix+"flag").c_str(), flag);
	// tree->SetBranchAddress((prefix+"merge_tag").c_str(), merge_tag);
	tree->SetBranchAddress((prefix+"x").c_str(), x);
	tree->SetBranchAddress((prefix+"y").c_str(), y);
	tree->SetBranchAddress((prefix+"z").c_str(), z);
	tree->SetBranchAddress((prefix+"energy").c_str(), energy);
	tree->SetBranchAddress((prefix+"time").c_str(), time);
	tree->SetBranchAddress((prefix+"time_flag").c_str(), time_flag);
}


void DssdMergeEvent::SetupOutput(TTree *tree) {
	tree->Branch("hit", &hit, "hit/s");
	tree->Branch("case", &case_tag, "case/i");
	tree->Branch("flag", flag, "flag[hit]/i");
	tree->Branch("merge_tag", merge_tag, "mtag[hit]/s");
	tree->Branch("x", x, "x[hit]/D");
	tree->Branch("y", y, "y[hit]/D");
	tree->Branch("z", z, "z[hit]/D");
	tree->Branch("energy", energy, "e[hit]/D");
	tree->Branch("time", time, "t[hit]/D");
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
	tree->SetBranchAddress((prefix+"time").c_str(), time);
	// tree->SetBranchAddress((prefix+"decode_entry").c_str(), decode_entry);
}


void AdssdMergeEvent::SetupOutput(TTree *tree) {
	tree->Branch("hit", &hit, "hit/s");
	tree->Branch("radius", radius, "r[hit]/D");
	tree->Branch("theta", theta, "theta[hit]/D");
	tree->Branch("phi", phi, "phi[hit]/D");
	tree->Branch("energy", energy, "e[hit]/D");
	tree->Branch("time", time, "time[hit]/D");
	// tree->Branch("decode_entry", decode_entry, "de[hit]/L");
}


}	// namespace ribll