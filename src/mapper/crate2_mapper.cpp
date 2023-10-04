#include "include/mapper/crate2_mapper.h"

namespace ribll {

Crate2Mapper::Crate2Mapper(unsigned int run)
: XiaMapper(run) {
}


int Crate2Mapper::Map(bool threshold) {
	TString input_file_name;
	input_file_name.Form
	(
		"%s%s_R%04d.root",
		kCrate2Path,  kCrate2FileName, run_
	);
	TTree *ipt = Initialize(input_file_name.Data());
	if (!ipt) {
		std::cerr << "Error: initialize tree from " << input_file_name
			<< " failed.\n";
		return -1;
	}

	// create output trees
	size_t t0d1_index = CreateOutputTree("t0d1", threshold);

	// show process
	printf("Mapping crate 2   0%%");
	fflush(stdout);
	long long entries = ipt->GetEntries();
	long long entry100 = entries / 100;

	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);

		energy_ = raw_energy_;
		timestamp_ = CalculateTimestamp(rate_, ts_);
		time_ = CalculateTime(rate_, timestamp_, cfd_, cfds_, cfdft_);
		decode_entry_ = entry;

		side_ = sid_ < 6 ? 0 : 1;
		switch ((sid_ - 2) % 4) {
			case 0:
				strip_ = ch_ + 32;
				break;
			case 1:
				strip_ = ch_ +  48;
				break;
			case 2:
				strip_ = ch_;
				break;
			default:
				strip_ = ch_ + 16;
				break;
		}

		if (side_ == 1 && (strip_ < 16 || strip_ >= 32)) {
			if (!threshold || raw_energy_ > 200) {
				FillTree(t0d1_index);
			}
		} else {
			if (!threshold || raw_energy_ > 300) {
				FillTree(t0d1_index);
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");


	for (size_t i = 0; i < opfs_.size(); ++i) {
		opfs_[i]->cd();
		opts_[i]->Write();
		opfs_[i]->Close();
	}

	MapStatistics statistics(run_, 2, threshold);
	statistics.Write();
	statistics.Print();

	return 0;
}

}		// namespace ribll