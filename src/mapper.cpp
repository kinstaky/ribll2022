#include "include/mapper.h"

#include <iostream>

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
	file_name.Form("%s%s%s-residual-%04d.root", kGenerateDataPath, kDecodeDir, name, run_);
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



//-----------------------------------------------------------------------------
//								crate 1 mapper
//-----------------------------------------------------------------------------

 Crate1Mapper::Crate1Mapper(int run)
 : XiaMapper(run) {
 }


int Crate1Mapper::Mapping() {
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
	size_t vme_trigger_index = CreateOutputTree("vt");
	size_t xia_trigger_index = CreateOutputTree("xt");
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

int Crate2Mapper::Mapping() {
	

	TString input_file_name;
	input_file_name.Form("%s%s_R%04d.root", kCrate2Path,  kCrate2FileName, run_);
	TTree *ipt = Initialize(input_file_name.Data());
	if (!ipt) {
		std::cerr << "Error: initialize tree from " << input_file_name << "failed.\n";
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
			FillTree(t0d2_index);
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