#include "include/mapper/crate3_mapper.h"

#include "include/statistics/map_statistics.h"

namespace ribll {

Crate3Mapper::Crate3Mapper(unsigned int run)
: run_(run) {
}


TTree* Crate3Mapper::Initialize(const char *file_name) {
	TFile *ipf = new TFile(file_name, "read");
	if (!ipf) {
		std::cerr << "Error: open file " << file_name << " failed.\n";
		return nullptr;
	}
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: get tree from " << file_name << " failed.\n";
		return nullptr;
	}

	ipt->SetBranchAddress("adc", adc_);
	ipt->SetBranchAddress("madc", madc_);
	ipt->SetBranchAddress("gdc", gdc_);
	ipt->SetBranchAddress("gmulti", gmulti_);

	return ipt;
}


size_t Crate3Mapper::CreatePPACTree(bool independent) {
	// output file name
	TString file_name;
	file_name.Form(
		"%s%svppac-%s-%04d.root",
		kGenerateDataPath,
		independent ? kFundamentalDir : kMappingDir,
		independent ? "fundamental" : "map",
		run_
	);
	opfs_.push_back(new TFile(file_name, "recreate"));

	TTree *opt = new TTree("tree", "tree of vppac");
	ppac_event_.SetupOutput(opt);
	if (!independent) {
		opt->Branch("timestamp", &align_time_, "ts/L");
	}

	opts_.push_back(opt);
	return opts_.size() - 1;
}


size_t Crate3Mapper::CreateADSSDTree(const char *name, bool independent) {
	TString file_name;
	file_name.Form(
		"%s%s%s-%s-%04d.root",
		kGenerateDataPath,
		independent ? kFundamentalDir : kMappingDir,
		name,
		independent ? "fundamental" : "map",
		run_
	);
	opfs_.push_back(new TFile(file_name, "recreate"));

	TTree *opt = new TTree("tree", TString::Format("tree of %s", name));
	dssd_event_.SetupOutput(opt);
	if (!independent) {
		opt->Branch("timestamp", &align_time_, "ts/L");
	}

	opts_.push_back(opt);
	return opts_.size() - 1;
}


size_t Crate3Mapper::CreateTofTree(bool independent) {
	TString file_name;
	file_name.Form(
		"%s%svtof-%s-%04d.root",
		kGenerateDataPath,
		independent ? kFundamentalDir : kMappingDir,
		independent ? "fundamental" : "map",
		run_
	);
	opfs_.push_back(new TFile(file_name, "recreate"));

	TTree *opt = new TTree("tree", "tree of vtof");
	tof_event_.SetupOutput(opt);
	if (!independent) {
		opt->Branch("timestamp", &align_time_, "ts/L");
	}

	opts_.push_back(opt);
	return opts_.size() - 1;
}


int Crate3Mapper::Map(bool independent) {
	// input decode file name
	TString input_file_name;
	input_file_name.Form(
		"%s%s%04d.root",
		kCrate3Path, kCrate3FileName, run_
	);
	// pointer to input file
	TTree *ipt = Initialize(input_file_name.Data());
	if (!ipt) {
		std::cerr << "Error: initialize tree from " << input_file_name
			<< " failed.\n";
		return -1;
	}

	// align file name
	TString align_file_name;
	align_file_name.Form(
		"%s%salign-%04d.root",
		kGenerateDataPath, kAlignDir, run_
	);
	// align-gdc file name
	TString align_gdc_file_name;
	align_gdc_file_name.Form(
		"%s%salign-gdc-%04u.root",
		kGenerateDataPath, kAlignDir, run_
	);

	if (!independent) {
		// add friend to combine with align tree
		if (!ipt->AddFriend("at=tree", align_file_name)) {
			std::cerr << "Error: Add friend from "
				<< align_file_name << " failed.\n";
			return -1;
		}
		// add friend to combine with align-gdc tree
		if (!ipt->AddFriend("ag=tree", align_gdc_file_name)) {
			std::cerr << "Error: Add friend from "
				<< align_gdc_file_name << "failed.\n";
			return -1;
		}
		// setup align tree branch
		ipt->SetBranchAddress("at.xia_time", &align_time_);
		// setup align-gdc tree branch
		ipt->SetBranchAddress("ag.gmulti", align_gmulti_);
		ipt->SetBranchAddress("ag.gdc", align_gdc_);
	}

	// create output trees
	size_t vtof_index = CreateTofTree(independent);
	size_t vppac_index = CreatePPACTree(independent);
	size_t tafd_index[2] = {
		CreateADSSDTree("tafd0", independent),
		CreateADSSDTree("tafd1", independent)
	};
	size_t tabd_index[6];
	for (size_t i = 0; i < 6; ++i) {
		tabd_index[i] = CreateADSSDTree(
			("tabd" + std::to_string(i)).c_str(),
			independent
		);
	}


	// total number of entries in input tree
	long long entries = ipt->GetEntries();
	// 1/100 of input entry, for showing process
	long long entry100 = entries / 100 + 1;
	// show start
	printf("Mapping crate 3   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		// read input event and align time
		ipt->GetEntry(entry);
		if (align_time_ < 0) continue;

		// tof
		tof_event_.time[0] = -1e5;
		tof_event_.time[1] = -1e5;
		tof_event_.cfd_flag = 0;
		for (size_t i = 0; i < 2; ++i) {
			if (gmulti_[0][125+i] > 0) {
				tof_event_.time[i] =
					(gdc_[0][125+i][0] - gdc_[0][127][0]) * 0.1;
			}
		}
		FillTree(vtof_index);

		// ppac
		ppac_event_.flag = 0;
		ppac_event_.cfd_flag = 0;
		ppac_event_.hit = ppac_event_.x_hit = ppac_event_.y_hit = 0;
		for (size_t i = 0; i < ppac_num; ++i) {
			ppac_event_.x1[i] = -1e5;
			ppac_event_.x2[i] = -1e5;
			ppac_event_.y1[i] = -1e5;
			ppac_event_.y2[i] = -1e5;
			ppac_event_.anode[i] = -1e5;

			// check existence of all signal of ppac
			for (size_t j = 0; j < 5; ++j) {
				if (gmulti_[0][i*5+j+96] > 0) {
					++ppac_event_.hit;
					ppac_event_.flag |= (1 << (i*5+j));
				}
			}
			// check x1 and x2
			if (gmulti_[0][i*5+96] > 0) {
				ppac_event_.x1[i] =
					(gdc_[0][i*5+96][0] - gdc_[0][127][0]) * 0.1;
			}
			if (gmulti_[0][i*5+97] > 0) {
				ppac_event_.x2[i] =
					(gdc_[0][i*5+97][0] - gdc_[0][127][0]) * 0.1;
			}
			if (gmulti_[0][i*5+96] > 0 && gmulti_[0][i*5+97] > 0) {
				++ppac_event_.x_hit;
			}

			// check y1 and y2
			if (gmulti_[0][i*5+98] > 0) {
				ppac_event_.y1[i] =
					(gdc_[0][i*5+98][0] - gdc_[0][127][0]) * 0.1;
			}
			if (gmulti_[0][i*5+99] > 0) {
				ppac_event_.y2[i] =
					(gdc_[0][i*5+99][0] - gdc_[0][127][0]) * 0.1;
			}
			if (gmulti_[0][i*5+98] > 0 && gmulti_[0][i*5+99] > 0) {
				++ppac_event_.y_hit;
			}

			// check anode
			if (gmulti_[0][i*5+100] > 0) {
				ppac_event_.anode[i] =
					(gdc_[0][i*5+100][0] - gdc_[0][127][0]) * 0.1;
			}
		}
		FillTree(vppac_index);

		// taf
		dssd_event_.cfd_flag = 0;
		for (size_t i = 0; i < 2; ++i) {
			dssd_event_.front_hit = dssd_event_.back_hit = 0;
			// check front side
			for (size_t j = 0; j < 16; ++j) {
				double energy
					= madc_[vtaf_front_module[i]][vtaf_front_channel[i]+j];
				if (align_gmulti_[i*16+j] > 0 && energy > 200) {
					size_t index = dssd_event_.front_hit;
					if (index >= 8) break;
					dssd_event_.front_strip[index] = j;
					dssd_event_.front_energy[index] = energy;
					if (!independent) {
						dssd_event_.front_time[index]
							= (align_gdc_[i*16+j][0] - align_gdc_[127][0]) * 0.1;
					} else {
						dssd_event_.front_time[index] = 0.0;
					}
					++dssd_event_.front_hit;
				}
			}
			// check back side
			for (size_t j = 0; j < 8; ++j) {
				double energy
					= madc_[vtaf_back_module[i]][vtaf_back_channel[i]+j];
				if (energy > 200) {
					size_t index = dssd_event_.back_hit;
					if (index >= 8) break;
					dssd_event_.back_strip[index] = j;
					dssd_event_.back_energy[index] = energy;
					++dssd_event_.back_hit;
				}
			}
			if (
				dssd_event_.front_hit > 0 && dssd_event_.front_hit <= 8
				&& dssd_event_.back_hit > 0 && dssd_event_.back_hit <= 8
			) {
				FillTree(tafd_index[i]);
			}
		}

		// tab
		for (size_t i = 0; i < 6; ++i) {
			dssd_event_.front_hit = dssd_event_.back_hit = 0;
			// check front side
			for (size_t j = 0; j < 16; ++j) {
				double energy
					= adc_[tab_front_module[i]][tab_front_channel[i]+j];
				if (gmulti_[0][((6-i)%6)*16+j] > 0 && energy > 200) {
					size_t index = dssd_event_.front_hit;
					if (index >= 8) break;
					dssd_event_.front_strip[index] = j;
					dssd_event_.front_energy[index] = energy;
					dssd_event_.front_time[index]
						= (gdc_[0][((6-i)%6)*16+j][0] - gdc_[0][127][0]) * 0.1;
					++dssd_event_.front_hit;
				}
			}
			// check back side
			for (size_t j = 0; j < 8; ++j) {
				double energy
					= adc_[tab_back_module[i]][tab_back_channel[i]+j];
				if (energy > 200) {
					size_t index = dssd_event_.back_hit;
					if (index >= 8) break;
					dssd_event_.back_strip[index] = j;
					dssd_event_.back_energy[index] = energy;
					++dssd_event_.back_hit;
				}

			}
			if (
				dssd_event_.front_hit > 0 && dssd_event_.front_hit <= 8
				&& dssd_event_.back_hit > 0 && dssd_event_.back_hit <= 8
			) {
				FillTree(tabd_index[i]);
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

	MapStatistics statistics(run_, 3, true);
	statistics.Write();
	statistics.Print();

	return 0;
}

} 		// namespace ribll