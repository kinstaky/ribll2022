#include "include/mapper.h"

#include <iostream>
#include <fstream>

#include <TString.h>

#include "include/defs.h"

namespace ribll {


	
//-----------------------------------------------------------------------------
//								xia mapper
//-----------------------------------------------------------------------------

XiaMapper::XiaMapper(int run)
: run_(run) {
}


XiaMapper::~XiaMapper() {
	for (TFile* ipf : ipfs_) {
		ipf->Close();
	}
}


TTree* XiaMapper::Initialize(const char *file_name) {
	TFile *input_file = new TFile(file_name, "read");
	if (!input_file) {
		std::cerr << "Error: open file " << file_name << " failed.\n";
		return nullptr;
	}
	ipfs_.push_back(input_file);

	TTree *ipt = (TTree*)input_file->Get("tree");
	if (!ipt) {
		std::cerr << "Error: get tree from " << file_name << " failed.\n";
		return nullptr;
	}

	ipt->SetBranchAddress("sr", &rate_);
	ipt->SetBranchAddress("sid", &sid_);
	ipt->SetBranchAddress("ch", &ch_);
	ipt->SetBranchAddress("ts", &ts_);
	ipt->SetBranchAddress("cfd", &cfd_);
	ipt->SetBranchAddress("cfds", &cfds_);
	ipt->SetBranchAddress("cfdft", &cfdft_);
	ipt->SetBranchAddress("evte", &raw_energy_);
	return ipt;
}


Long64_t XiaMapper::CalculateTimestamp(Short_t rate, Long64_t ts) {
	Long64_t result = 0;
	switch (rate) {
		case 100:
		case 500:
			result = ts * 10;
			break;
		case 250:
			result = ts * 8;
			break;
		default:
			;
	}
	return result;
}


Double_t XiaMapper::CalculateTime(
	Short_t rate,
	Long64_t timestamp,
	Short_t cfd,
	Short_t cfds,
	Bool_t cfdft
) {
	Double_t result = 0.0;
	switch (rate) {
		case 100:
			result = timestamp
				+ Double_t(cfd) / 32768.0 * 10.0;
			break;			
		case 250:
			result = timestamp
				+ (cfdft ? 0.0 : (Double_t(cfd) / 166384.0 - cfds) * 4.0);
			break;
		case 500:
			result = timestamp
				+ (cfdft ? 0.0 : (Double_t(cfd) / 8192.0 + cfds - 1) * 2.0);
			break;
		default:
			result = 0.0;
	}
	return result;
}


size_t XiaMapper::CreateOutputTree(const char *name) {
	TString file_name;
	file_name.Form("%s%s%s-map-%04d.root", kGenerateDataPath, kMappingDir, name, run_);
	opfs_.push_back(new TFile(file_name, "recreate"));

	TTree *opt = new TTree("tree", TString::Format("tree of %s", name));
	opt->Branch("index", &detector_index_, "index/s");
	opt->Branch("side", &side_, "side/s");
	opt->Branch("strip", &strip_, "s/s");
	opt->Branch("timestamp", &timestamp_, "ts/L");
	opt->Branch("time", &time_, "t/D");
	opt->Branch("energy", &energy_, "e/D");

	opts_.push_back(opt);
	return opts_.size() - 1;
}

size_t XiaMapper::CreateResidualTree(const char *name) {
	TString file_name;
	file_name.Form("%s%s%s-residual-%04d.root", kGenerateDataPath, kMappingDir, name, run_);
	opfs_.push_back(new TFile(file_name, "recreate"));

	TTree *opt = new TTree("tree", TString::Format("residual tree of %s", name));
	opt->Branch("sr", &rate_);
	opt->Branch("sid", &sid_);
	opt->Branch("ch", &ch_);
	opt->Branch("ts", &ts_);
	opt->Branch("cfd", &cfd_);
	opt->Branch("cfds", &cfds_);
	opt->Branch("cfdft", &cfdft_);
	opt->Branch("evte", &raw_energy_);

	opts_.push_back(opt);
	return opts_.size() - 1;
}


size_t XiaMapper::CreateTriggerTree(const char *name) {
	TString file_name;
	file_name.Form("%s%s%s-map-%04d.root", kGenerateDataPath, kMappingDir, name, run_);
	opfs_.push_back(new TFile(file_name, "recreate"));

	TTree *opt = new TTree("tree", TString::Format("tree of %s", name));
	opt->Branch("timestamp", &timestamp_, "ts/L");

	opts_.push_back(opt);
	return opts_.size() - 1;
}


//-----------------------------------------------------------------------------
//								crate 1 mapper
//-----------------------------------------------------------------------------

 Crate1Mapper::Crate1Mapper(int run)
 : XiaMapper(run) {
 }


int Crate1Mapper::Map() {
	TString input_file_name;
	input_file_name.Form("%s%s_R%04d.root", kCrate1Path, kCrate1FileName, run_);
	TTree *ipt = Initialize(input_file_name.Data());
	if (!ipt) {
		std::cerr << "Error: intialize tree from " << input_file_name << " failed\n";
		return -1;
	}

	// create residual tree
	size_t residual_index = CreateResidualTree("c1");

	// create output trees
	size_t tof_index = CreateOutputTree("tof");
	size_t vme_trigger_index = CreateTriggerTree("vt");
	size_t xia_trigger_index = CreateTriggerTree("xt");
	size_t xia_ppac_index = CreateOutputTree("xppac");
	size_t xtaf_index = CreateOutputTree("xtaf");
	size_t t0ssd_index = CreateOutputTree("t0ssd");
	size_t t1ssd_index = CreateOutputTree("t1ssd");
	size_t t0csi_index = CreateOutputTree("t0csi");
	size_t t1csi_index = CreateOutputTree("t1csi");
	size_t tafcsi_index = CreateOutputTree("tafcsi");
	size_t tabcsi_index = CreateOutputTree("tabcsi");

	// show process
	printf("Mapping crate 1   0%%");
	fflush(stdout);
	Long64_t entry100 = ipt->GetEntries() / 100;
	
	for (Long64_t entry = 0; entry < 10'000'000; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);
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
					FillTree(t0ssd_index);
					break;
				case 4:
					// t0ssd 2, fall through
				case 5:
					// t0ssd 3
					detector_index_ = ch_ - 3;
					FillTree(t0ssd_index);
					break;
				case 6:
					// t1ssd
					FillTree(t1ssd_index);
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
					// tabcsi-0A
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
				detector_index_ = ch_;
				FillTree(tabcsi_index);
			} else {
				// t0 csi
				detector_index_ = ch_ - 12;
				FillTree(t0csi_index);
			}
		} else if (sid_ == 5) {
			if (ch_ < 12) {
				// taf csi
				detector_index_ = ch_;
				FillTree(tafcsi_index);
			} else {
				// t1 csi
				detector_index_ = ch_ - 12;
				FillTree(t1csi_index);
			}
		} else if (sid_ >= 6 && sid_ < 12) {
			// taf in xia
			switch (sid_) {
				case 6:
					detector_index_ = 4;
					strip_ = ch_;
					break;
				case 7:
					detector_index_ = 5;
					strip_ = ch_;
					break;
				case 8:
					detector_index_ = ch_ < 8 ? 4 : 5;
					side_ = 1;
					strip_ = 7 - (ch_ % 8);
					break;
				case 9:
					detector_index_ = 2;
					strip_ = ch_;
					break;
				case 10:
					detector_index_ = 3;
					strip_ = ch_;
					break;
				case 11:
					detector_index_ = ch_ < 8 ? 2 : 3;
					side_ = 1;
					strip_ = 7 - (ch_ % 8);
					break;
				default:
					std::cerr << "Error: switch taf should not be here.\n";
					return -1;
			}
			FillTree(xtaf_index);
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
	return 0;
}



//-----------------------------------------------------------------------------
//								crate 2 mapper
//-----------------------------------------------------------------------------

Crate2Mapper::Crate2Mapper(int run)
:XiaMapper(run) {
}

int Crate2Mapper::Map() {
	
	TString input_file_name;
	input_file_name.Form("%s%s_R%04d.root", kCrate2Path,  kCrate2FileName, run_);
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
	printf("Mapping crate 2   0%%");
	fflush(stdout);
	Long64_t entry100 = ipt->GetEntries() / 100;
	for (Long64_t entry = 0; entry < ipt->GetEntries(); ++entry) {
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
			FillTree(t1d1_index);
		} else if (sid_ < 10) {
			// t0d3
			detector_index_ = 0;
			side_ = sid_ < 8 ? 0 : 1;
			strip_ = sid_ % 2 == 0 ? ch_ : ch_ + 16;
			FillTree(t0d3_index);
		} else {
			// t0d2
			detector_index_ = 0;
			side_ = sid_ < 12 ? 0 : 1;
			strip_ = sid_ % 2 == 0 ? ch_ : ch_ + 16;
			if (raw_energy_ > 20) {
				FillTree(t0d2_index);
			}
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");



	// read events of t0d3 recorded in crate1
	input_file_name.Form("%s%sc1-residual-%04d.root", kGenerateDataPath, kMappingDir, run_);
	TTree *res_ipt = Initialize(input_file_name.Data());
	if (!res_ipt) {
		std::cerr << "Error: initialize tree from " << input_file_name << " failed.\n";
		return -1;
	}

	// show process
	printf("Mapping crate 1 residual   0%%");
	fflush(stdout);
	entry100 = res_ipt->GetEntries() / 100;
	for (Long64_t entry = 0; entry < res_ipt->GetEntries(); ++entry) {
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
		FillTree(t0d3_index);
	}
	// show finish
	printf("\b\b\b\b100%%\n");


	for (size_t i = 0; i < opfs_.size(); ++i) {
		opfs_[i]->cd();
		opts_[i]->Write();
		opfs_[i]->Close();
	}
	return 0;
}


//-----------------------------------------------------------------------------
//								crate 3 mapper
//-----------------------------------------------------------------------------

Crate3Mapper::Crate3Mapper(int run)
: XiaMapper(run) {
}


int Crate3Mapper::Map() {
	TString input_file_name;
	input_file_name.Form("%s%s_R%04d.root", kCrate3Path,  kCrate3FileName, run_);
	TTree *ipt = Initialize(input_file_name.Data());
	if (!ipt) {
		std::cerr << "Error: initialize tree from " << input_file_name << " failed.\n";
		return -1;
	}

	// create output trees
	size_t t0d1_index = CreateOutputTree("t0d1");

	// show process
	printf("Mapping crate 3   0%%");
	fflush(stdout);
	Long64_t entry100 = ipt->GetEntries() / 100;
	
	for (Long64_t entry = 0; entry < ipt->GetEntries(); ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);
		energy_ = raw_energy_;
		timestamp_ = CalculateTimestamp(rate_, ts_);
		time_ = CalculateTime(rate_, timestamp_, cfd_, cfds_, cfdft_);
		
		side_ = sid_ < 6 ? 1 : 0;
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

		if (raw_energy_ > 20) {
			FillTree(t0d1_index);
		}
	}
	// show finish
	printf("\b\b\b\b100%%\n");


	for (size_t i = 0; i < opfs_.size(); ++i) {
		opfs_[i]->cd();
		opts_[i]->Write();
		opfs_[i]->Close();
	}
	return 0;
}


//-----------------------------------------------------------------------------
//								crate 4 mapper
//-----------------------------------------------------------------------------

Crate4Mapper::Crate4Mapper(int run)
: run_(run) {
}


TTree* Crate4Mapper::Initialize(const char *file_name) {
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


size_t Crate4Mapper::CreatePPACTree(const char *name) {
	TString file_name;
	file_name.Form("%s%s%s-corr-%04d.root", kGenerateDataPath, kCorrelationDir, name, run_);
	opfs_.push_back(new TFile(file_name, "recreate"));

	TTree *opt = new TTree("tree", TString::Format("tree of %s", name));
	opt->Branch("timestamp", &ppac_event_.timestamp, "ts/L");
	opt->Branch("flag", &ppac_event_.flag, "flag/I");
	opt->Branch("hit", &ppac_event_.hit, "hit/s");
	opt->Branch("xhit", &ppac_event_.x_hit, "xhit/s");
	opt->Branch("yhit", &ppac_event_.y_hit, "yhit/s");
	opt->Branch("x", ppac_event_.x, TString::Format("x[%llu]/D", ppac_num).Data());
	opt->Branch("y", ppac_event_.y, TString::Format("y[%llu]/D", ppac_num).Data());

	opts_.push_back(opt);
	return opts_.size() - 1;
}


size_t Crate4Mapper::CreateADSSDTree(const char *name) {
	TString file_name;
	file_name.Form("%s%s%s-corr-%04d.root", kGenerateDataPath, kCorrelationDir, name, run_);
	opfs_.push_back(new TFile(file_name, "recreate"));

	TTree *opt = new TTree("tree", TString::Format("tree of %s", name));
	opt->Branch("timestamp", &dssd_event_.timestamp, "ts/L");
	opt->Branch("index", &dssd_event_.index, "index/s");
	opt->Branch("xhit", &dssd_event_.x_hit, "xhit/s");
	opt->Branch("yhit", &dssd_event_.y_hit, "yhit/s");
	opt->Branch("x_strip", dssd_event_.x_strip, "xs[xhit]/s");
	opt->Branch("y_strip", dssd_event_.y_strip, "ys[yhit]/s");
	opt->Branch("x_time", dssd_event_.x_time, "xt[xhit]/D");
	opt->Branch("y_time", dssd_event_.y_time, "yt[yhit]/D");
	opt->Branch("x_energy", dssd_event_.x_energy, "xe[xhit]/D");
	opt->Branch("y_energy", dssd_event_.y_energy, "ye[yhit]/D");

	opts_.push_back(opt);
	return opts_.size() - 1;
}


int Crate4Mapper::Map() {
	TString input_file_name;
	input_file_name.Form("%s%s%04d.root", kCrate4Path, kCrate4FileName, run_);
	TTree *ipt = Initialize(input_file_name.Data());
	if (!ipt) {
		std::cerr << "Error: initialize tree from " << input_file_name << " failed.\n";
		return -1;
	}

	// create output trees
	size_t vppac_index = CreatePPACTree("vppac");
	size_t vtaf_index = CreateADSSDTree("vtaf");
	size_t tab_index = CreateADSSDTree("tab");

	// read align timestamp from file
	TString align_file_name;
	align_file_name.Form("%s%stimestamp-%04d.txt", kGenerateDataPath, kAlignDir, run_);
	std::ifstream align_in(align_file_name.Data());
	if (!align_in.good()) {
		std::cerr << "Error: open file " << align_file_name << " failed.\n";
		return -1;
	}

	// show process
	printf("Mapping crate 4   0%%");
	fflush(stdout);
	Long64_t entry100 = ipt->GetEntries() / 100;
	
	for (Long64_t entry = 0; entry < ipt->GetEntries(); ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);

		// read align timestamp
		Long64_t timestamp;
		align_in >> timestamp;
		if (timestamp < 0) {
			continue;
		}
		
		
		// ppac
		ppac_event_.timestamp = timestamp;
		ppac_event_.flag = 0;
		ppac_event_.hit = ppac_event_.x_hit = ppac_event_.y_hit = 0;
		for (size_t i = 0; i < ppac_num; ++i) {
			// check existence of all signal of ppac
			for (size_t j = 0; j < 5; ++j) {
				if (gmulti_[0][i*5+j+96] > 0) {
					++ppac_event_.hit;
					ppac_event_.flag |= (1 << (i*5+j)); 
				}
			}
			// check x1 and x2
			if (gmulti_[0][i*5+96] > 0 && gmulti_[0][i*5+97] > 0) {
				++ppac_event_.x_hit;
				ppac_event_.x[i] = double(gdc_[0][i*5+96][0] - gdc_[0][i*5+97][0]) * 0.1;
			} else {
				ppac_event_.x[i] = 0.0;
			}
			// check y1 and y2
			if (gmulti_[0][i*5+98] > 0 && gmulti_[0][i*5+99] > 0) {
				++ppac_event_.y_hit;
				ppac_event_.y[i] = double(gdc_[0][i*5+98][0] - gdc_[0][i*5+99][0]) * 0.1; 
			} else {
				ppac_event_.y[i] = 0.0;
			}
		}
		FillTree(vppac_index);

		// vtaf
		dssd_event_.timestamp = timestamp;
		for (size_t i = 0; i < 2; ++i) {
			dssd_event_.index = i;
			dssd_event_.x_hit = dssd_event_.y_hit = 0;
			// check front side
			for (size_t j = 0; j < 16; ++j) {
				double energy = madc_[vtaf_front_module[i]][vtaf_front_channel[i]+j];
				// if (energy > 30 && gmulti_[1][i*16+j] == 1) {
				if (energy > 200) {
					size_t index = dssd_event_.x_hit;
					dssd_event_.x_strip[index] = j;
					dssd_event_.x_energy[index] = energy;
					dssd_event_.x_time[index] = (gdc_[1][i*16+j][0] - gdc_[1][127][0]) * 0.1;
					++dssd_event_.x_hit;
				}
			}
			// check back side
			for (size_t j = 0; j < 8; ++j) {
				double energy = madc_[vtaf_back_module[i]][vtaf_back_channel[i]-j];
				if (energy > 200) {
					size_t index = dssd_event_.y_hit;
					dssd_event_.y_strip[index] = j;
					dssd_event_.y_energy[index] = energy;
					++dssd_event_.y_hit;
				}
			}
			if (dssd_event_.x_hit && dssd_event_.y_hit) {
				FillTree(vtaf_index);
			}
		}

		// tab
		for (size_t i = 0; i < 6; ++i) {
			dssd_event_.index = i;
			dssd_event_.x_hit = dssd_event_.y_hit = 0;
			// check front side
			for (size_t j = 0; j < 16; ++j) {
				double energy = adc_[tab_front_module[i]][tab_front_channel[i]+j];
				// if (energy > 30 && gmulti_[0][i*16+j] == 1) {
				if (energy > 200) {
					size_t index = dssd_event_.x_hit;
					dssd_event_.x_strip[index] = j;
					dssd_event_.x_energy[index] = energy;
					dssd_event_.x_time[index] = (gdc_[0][i*16+j][0] - gdc_[0][127][0]) * 0.1;
					++dssd_event_.x_hit;
				}
			}
			// check back side
			for (size_t j = 0; j < 8; ++j) {
				double energy = adc_[tab_back_module[i]][tab_back_channel[i]-j];
				if (i == 5) {
					// special case
					energy = adc_[4][j];
				}
				if (energy > 200) {
					size_t index = dssd_event_.y_hit;
					dssd_event_.y_strip[index] = j;
					dssd_event_.y_energy[index] = energy;
					++dssd_event_.y_hit;
				}
			}
			if (dssd_event_.x_hit && dssd_event_.y_hit) {
				FillTree(tab_index);
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
	return 0;
}

}