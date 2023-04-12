#include "include/detector/adssd.h"

#include <iostream>

namespace ribll {

Adssd::Adssd(
	unsigned int run,
	const std::string &name,
	const std::string &tag
)
: Dssd(run, name, tag) {
}


ROOT::Math::Polar3DVector Adssd::CalculatePosition(
	unsigned short front_strip,
	unsigned short back_strip
) const {
	double radius = (radius_range_.second-radius_range_.first) / FrontStrip();
	radius  = radius * (front_strip + 0.5) + radius_range_.first;
	double phi = (phi_range_.second - phi_range_.first) / BackStrip();
	phi = phi * (back_strip + 0.5) + phi_range_.first;
	ROOT::Math::Polar3DVector result(radius, 90.0*TMath::DegToRad(), phi);
	result += center_;
	return result;
}


int Adssd::Merge(double energy_diff) {
	// input file name
	TString fundamental_file_name;
	fundamental_file_name.Form(
		"%s%s%s-fundamental-%s%04u.root",
		kGenerateDataPath,
		kFundamentalDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// input file
	TFile *ipf = new TFile(fundamental_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< fundamental_file_name << " failed.\n";
	}
	// input event
	DssdFundamentalEvent fundamental_event;
	// setup input branches
	fundamental_event.SetupInput(ipt);
	// for convenient
	unsigned short &fhit = fundamental_event.front_hit;
	unsigned short &bhit = fundamental_event.back_hit;
	unsigned short *fs = fundamental_event.front_strip;
	unsigned short *bs = fundamental_event.back_strip;
	double *fe = fundamental_event.front_energy;
	double *be = fundamental_event.back_energy;


	// output file name
	TString merge_file_name;
	merge_file_name.Form(
		"%s%s%s-merge-%s%04u.root",
		kGenerateDataPath,
		kMergeDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile *opf = new TFile(merge_file_name, "recreate");
	// relatetive difference of front and back side energy
	TH1F *hrd = new TH1F(
		"hrd", "relateive difference of front and back side energy",
		1000, 0.0, 1.0
	);
	// output tree
	TTree *opt = new TTree("tree", "tree of merged events");
	// output event
	AdssdMergeEvent merge_event;
	// setup output branches
	merge_event.SetupOutput(opt);

	// read normalized parameters
	if (ReadNormalizeParameters()) {
		std::cerr << "Error: Read normalize parameters failed.\n";
		return -1;
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100;
	// show start
	printf("Writing merged events   0%%");
	fflush(stdout);
	// loop over events
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);
		// initialize
		merge_event.hit = 0;

		if (fhit == 1 && bhit == 1) {
			if (fe[0] < 1e4 && be[0] < 1e4) {
				fe[0] = NormEnergy(0, fs[0], fe[0]);
				be[0] = NormEnergy(1, bs[1], be[0]);
				double diff = RelativeDifference(fe[0], be[0]);
				hrd->Fill(diff);
				if (diff < energy_diff) {
					merge_event.energy[0] = fe[0];
					auto position = CalculatePosition(fs[0], bs[0]);
					merge_event.radius[0] = position.R();
					merge_event.phi[0] = position.Phi();
					merge_event.theta[0] = position.Theta();
					merge_event.hit = 1;
				}
			}
		}

		opt->Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save and close files
	hrd->Write();
	opt->Write();
	opf->Close();
	ipf->Close();

	return 0;
}

}		// namespace ribll