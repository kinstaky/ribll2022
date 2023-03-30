#include "include/mapper/crate0_mapper.h"

namespace ribll {

Crate0Mapper::Crate0Mapper(unsigned int run)
: XiaMapper(run) {
}


int Crate0Mapper::Map(bool threshold) {
	TString input_file_name;
	input_file_name.Form("%s%s_R%04d.root", kCrate0Path, kCrate0FileName, run_);
	TTree *ipt = Initialize(input_file_name.Data());
	if (!ipt) {
		std::cerr << "Error: intialize tree from " << input_file_name
			<< " failed\n";
		return -1;
	}

	// create residual tree
	size_t residual_index = CreateResidualTree("c0");

	// create output trees
	size_t tof_index = CreateOutputTree("tof", threshold);
	size_t vme_trigger_index = CreateTriggerTree("vt", threshold);
	size_t xia_trigger_index = CreateTriggerTree("xt", threshold);
	size_t xia_ppac_index = CreateOutputTree("xppac", threshold);
	size_t taf2_index = CreateOutputTree("taf2", threshold);
	size_t taf3_index = CreateOutputTree("taf3", threshold);
	size_t taf4_index = CreateOutputTree("taf4", threshold);
	size_t taf5_index = CreateOutputTree("taf5", threshold);
	size_t t0s1_index = CreateOutputTree("t0s1", threshold);
	size_t t0s2_index = CreateOutputTree("t0s2", threshold);
	size_t t0s3_index = CreateOutputTree("t0s3", threshold);
	size_t t1s1_index = CreateOutputTree("t1s1", threshold);
	size_t t0csi_index = CreateOutputTree("t0csi", threshold);
	size_t t1csi_index = CreateOutputTree("t1csi", threshold);
	size_t tafcsi_index = CreateOutputTree("tafcsi", threshold);
	size_t tabcsi_index = CreateOutputTree("tabcsi", threshold);

	const unsigned short tabcsi_channel_index[12] = {
		0, 1, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3
	};

	// show process
	printf("Mapping crate 0   0%%");
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
		if (raw_energy_ == 0) continue;

		energy_ = raw_energy_;
		timestamp_ = CalculateTimestamp(rate_, ts_);
		time_ = CalculateTime(rate_, timestamp_, cfd_, cfds_, cfdft_);
		detector_index_ = 0;
		side_ = 0;
		strip_ = 0;

		if (sid_ == 2) {
			switch (ch_) {
				case 0:
					// tof 1, fall through
				case 1:
					// tof 2
					detector_index_ = ch_;
					FillTree(tof_index);
					break;
				case 2:
					// t0ssd 1
					FillTree(t0s1_index);
					break;
				case 4:
					// t0ssd 2
					FillTree(t0s2_index);
					break;
				case 5:
					// t0ssd 3
					FillTree(t0s3_index);
					break;
				case 6:
					// t1ssd
					FillTree(t1s1_index);
					break;
				case 8:
					// xia trigger
					FillTree(xia_trigger_index);
					break;
				case 9:
					// vme trigger
					FillTree(vme_trigger_index);
					break;
				case 15:
					// tabcsi-0
					FillTree(tabcsi_index);
					break;
				default:
					FillTree(residual_index);
			}
		} else if (sid_ == 3) {
			// ppac in xia
			// set ppac index
			detector_index_ = ch_ > 0 ? (ch_ - 1) / 5 : 0;
			// set ppac side and strip, side 0 for X1, side 1 for X2,
			// side 2 for Y1, side 3 for Y2, side 4 for A
			side_ = ch_ > 0 ? (ch_ - 1) % 5 : 0;
			FillTree(xia_ppac_index);
		} else if (sid_ == 4) {
			if (ch_ == 0) {
				continue;
			} else if (ch_ < 12) {
				// tab csi
				detector_index_ = tabcsi_channel_index[ch_];
				FillTree(tabcsi_index);
			} else {
				// t0 csi
				detector_index_ = ch_ - 12;
				FillTree(t0csi_index);
			}
		} else if (sid_ == 5) {
			if (ch_ < 12) {
				switch (ch_) {
					// special case for 8 and 9
					case 8:
						detector_index_ = 9;
						break;
					case 9:
						detector_index_ = 8;
						break;
					default:
						detector_index_ = ch_;
				}
				// taf csi
				FillTree(tafcsi_index);
			} else {
				// t1 csi
				detector_index_ = ch_ - 12;
				FillTree(t1csi_index);
			}
		} else if (sid_ == 6) {
			// taf4 front side
			strip_ = ch_;
			if (!threshold || raw_energy_ > 200) {
				FillTree(taf4_index);
			}
		} else if (sid_ == 7) {
			// taf5 front side
			strip_ = ch_;
			if (raw_energy_ > 100) {
				FillTree(taf5_index);
			}
		} else if (sid_ == 8) {
			// taf4/taf5 back side
			side_ = 1;
			if (ch_ < 8) {
				// taf5 back side
				strip_ = ch_;
				switch (strip_) {
					case 7:
						if (raw_energy_ > 600) FillTree(taf5_index);
						break;
					case 6:
						if (raw_energy_ > 500) FillTree(taf5_index);
						break;
					case 0:
					case 1:
					case 5:
						if (raw_energy_ > 450) FillTree(taf5_index);
						break;
					default:
						if (raw_energy_ > 400) FillTree(taf5_index);
				}
			} else {
				// taf4 back side
				strip_ = ch_ - 8;
				if (
					!threshold
					|| (
						(strip_ < 4 && raw_energy_ > 400)
						|| (strip_ >= 4 && strip_ < 6 && raw_energy_ > 1000)
						|| (strip_ >= 6 && raw_energy_ > 1200)
					)
				) {
					FillTree(taf4_index);
				}
			}
		} else if (sid_ == 9) {
			// taf2 front side
			strip_ = ch_;
			if (!threshold || raw_energy_ > 150) {
				FillTree(taf2_index);
			}
		} else if (sid_ == 10) {
			// taf3 front side
			strip_ = ch_;
			if (!threshold || raw_energy_ > 150) {
				FillTree(taf3_index);
			}
		} else if (sid_ == 11) {
			// taf2/taf3 back side
			side_ = 1;
			if (ch_ < 8) {
				// taf2 back side
				strip_ = ch_;
				if (!threshold || raw_energy_ > 200) {
					FillTree(taf2_index);
				}
			} else {
				// taf3 back side
				strip_ = ch_ - 8;
				if (!threshold || raw_energy_ > 400) {
					FillTree(taf3_index);
				}
			}
		} else {
			FillTree(residual_index);
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	for (size_t i = 0; i < opfs_.size(); ++i) {
		opfs_[i]->cd();
		opts_[i]->Write();
		opfs_[i]->Close();
	}

	MapStatistics statistics(run_, 0, threshold);
	statistics.Write();
	statistics.Print();

	return 0;
}

}		// namespace ribll