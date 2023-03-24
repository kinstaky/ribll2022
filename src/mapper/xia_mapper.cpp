#include "include/mapper/xia_mapper.h"

#include <iostream>

#include <TString.h>

#include "include/defs.h"

namespace ribll {



//-----------------------------------------------------------------------------
//								xia mapper
//-----------------------------------------------------------------------------

XiaMapper::XiaMapper(unsigned int run)
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


long long XiaMapper::CalculateTimestamp(short rate, long long ts) {
	long long result = 0;
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


double XiaMapper::CalculateTime(
	short rate,
	long long timestamp,
	short cfd,
	short cfds,
	bool cfdft
) {
	double result = 0.0;
	switch (rate)
	{
		case 100:
			result = timestamp
				+ double(cfd) / 32768.0 * 10.0;
			break;
		case 250:
			result = timestamp
				+ (cfdft ? 0.0 : (double(cfd) / 166384.0 - cfds) * 4.0);
			break;
		case 500:
			result = timestamp
				+ (cfdft ? 0.0 : (double(cfd) / 8192.0 + cfds - 1) * 2.0);
			break;
		default:
			result = 0.0;
	}
	return result;
}


size_t XiaMapper::CreateOutputTree(const char *name) {
	TString file_name;
	file_name.Form(
		"%s%s%s-map-%04d.root",
		kGenerateDataPath, kMappingDir, name, run_
	);
	opfs_.push_back(new TFile(file_name, "recreate"));

	TTree *opt = new TTree("tree", TString::Format("tree of %s", name));
	opt->Branch("index", &detector_index_, "index/s");
	opt->Branch("side", &side_, "side/s");
	opt->Branch("strip", &strip_, "s/s");
	opt->Branch("timestamp", &timestamp_, "ts/L");
	opt->Branch("time", &time_, "t/D");
	opt->Branch("energy", &energy_, "e/D");
	opt->Branch("cfd", &cfdft_, "cfd/O");

	opts_.push_back(opt);
	return opts_.size() - 1;
}

size_t XiaMapper::CreateResidualTree(const char *name) {
	TString file_name;
	file_name.Form(
		"%s%s%s-residual-%04d.root",
		kGenerateDataPath, kMappingDir, name, run_
	);
	opfs_.push_back(new TFile(file_name, "recreate"));

	TTree *opt = new TTree(
		"tree",
		TString::Format("residual tree of %s", name)
	);
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
	file_name.Form(
		"%s%s%s-map-%04d.root",
		kGenerateDataPath, kMappingDir, name, run_
	);
	opfs_.push_back(new TFile(file_name, "recreate"));

	TTree *opt = new TTree("tree", TString::Format("tree of %s", name));
	opt->Branch("time", &time_, "time/D");
	opt->Branch("timestamp", &timestamp_, "ts/L");
	opt->Branch("cfd", &cfdft_, "cfd/O");

	opts_.push_back(opt);
	return opts_.size() - 1;
}

} 		// namespace ribll