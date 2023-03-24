#include "include/mapper/crate1_mapper.h"

namespace ribll {

Crate1Mapper::Crate1Mapper(unsigned int run)
:XiaMapper(run) {
}

int Crate1Mapper::Map(bool threshold) {
	TString input_file_name;
	input_file_name.Form(
		"%s%s_R%04d.root",
		kCrate1Path,  kCrate1FileName, run_
	);
	TTree *ipt = Initialize(input_file_name.Data());
	if (!ipt) {
		std::cerr << "Error: initialize tree from " << input_file_name << " failed.\n";
		return -1;
	}

	// create output trees
	size_t t1d1_index = CreateOutputTree("t1d1");
	size_t t0d3_index = CreateOutputTree("t0d3");
	size_t t0d2_index = CreateOutputTree("t0d2");

	// show process
	printf("Mapping crate 1   0%%");
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

		if (sid_ < 6) {
			// t1d1
			detector_index_ = 0;
			side_ = sid_ < 4 ? 0 : 1;
			switch (sid_) {
				case 2:
					strip_ = 31 - ch_;
					break;
				case 3:
					strip_ = 15 - ch_;
					break;
				case 4:
					strip_ = ch_;
					break;
				case 5:
					strip_ = ch_ + 16;
					break;
				default:
					std::cerr << "Error: switch t1d1 should not be here.\n";
					return -1;
			}
			if (!threshold || raw_energy_ > 4000) FillTree(t1d1_index);
		}
		else if (sid_ < 10) {
			// t0d3
			detector_index_ = 0;
			side_ = sid_ < 8 ? 0 : 1;
			strip_ = sid_ % 2 == 0 ? ch_ : ch_ + 16;
			if (!side_) {
				if (!threshold || raw_energy_ > 100) {
					FillTree(t0d3_index);
				}
			} else {
				if (!threshold || raw_energy_ > 50) {
					FillTree(t0d3_index);
				}
			}
		} else {
			// t0d2
			detector_index_ = 0;
			side_ = sid_ < 12 ? 0 : 1;
			strip_ = sid_ % 2 == 0 ? ch_ : ch_ + 16;
			if (!side_) {
				if (!threshold || raw_energy_ > 300) {
					FillTree(t0d2_index);
				}
			} else {
				if (!threshold || raw_energy_ > 800) {
					FillTree(t0d2_index);
				}
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");


	// read events of t0d3 recorded in crate1
	input_file_name.Form(
		"%s%sc0-residual-%04d.root",
		kGenerateDataPath, kMappingDir, run_
	);
	TTree *res_ipt = Initialize(input_file_name.Data());
	if (!res_ipt) {
		std::cerr << "Error: initialize tree from " << input_file_name
			<< " failed.\n";
		return -1;
	}

	// show process
	printf("Mapping crate 0 residual   0%%");
	fflush(stdout);
	entries = res_ipt->GetEntries();
	entry100 = entries / 100;
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		res_ipt->GetEntry(entry);
		if (sid_ != 2) continue;


		energy_ = raw_energy_;
		timestamp_ = CalculateTimestamp(rate_, ts_);
		time_ = CalculateTime(rate_, timestamp_, cfd_, cfds_, cfdft_);

		switch (ch_) {
			case 10:
				side_ = 0;
				strip_ = 15;
				break;
			case 11:
				side_ = 1;
				strip_ = 15;
				break;
			case 14:
				side_ = 1;
				strip_ = 30;
				break;
			default:
				continue;
		}

		if (!side_) {
			if (!threshold || raw_energy_ > 100) {
				FillTree(t0d3_index);
			}
		} else {
			if (!threshold || raw_energy_ > 50) {
				FillTree(t0d3_index);
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

	MapStatistics statistics(run_, 1, threshold);
	statistics.Write();
	statistics.Print();

	return 0;
}

}		// namespace ribll