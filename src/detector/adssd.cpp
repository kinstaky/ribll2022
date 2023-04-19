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


int Adssd::Merge(double) {
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
	TFile ipf(fundamental_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf.Get("tree");
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
	TFile opf(merge_file_name, "recreate");
	// output tree
	TTree opt("tree", "tree of merged events");
	// output event
	AdssdMergeEvent merge_event;
	// setup output branches
	merge_event.SetupOutput(&opt);

	// calibrated parameters
	double cali_params[16][2];
	// calibration parameters file name
	TString param_file_name;
	param_file_name.Form(
		"%s%s%s-alpha-cali-param.txt",
		kGenerateDataPath,
		kCalibrationDir,
		name_.c_str()
	);
	// parameters file
	std::ifstream fin(param_file_name);
	if (!fin.good()) {
		std::cerr << "Error: Open file "
			<< param_file_name << " failed.\n";
		return -1;
	}
	// read parameters from file
	for (size_t i = 0; i < 16; ++i) {
		fin >> cali_params[i][0] >> cali_params[i][1];
		double tmp;
		fin >> tmp >> tmp >> tmp;
	}
	// close parameters file
	fin.close();


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
			if (fe[0] < 1e4) {
				fe[0] = cali_params[fs[0]][0] + cali_params[fs[0]][1] * fe[0];
				merge_event.energy[0] = fe[0];
				auto position = CalculatePosition(fs[0], bs[0]);
				merge_event.radius[0] = position.R();
				merge_event.phi[0] = position.Phi();
				merge_event.theta[0] = position.Theta();
				merge_event.hit = 1;
			}
		}

		opt.Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save and close files
	opt.Write();
	opf.Close();
	ipf.Close();

	return 0;
}

}		// namespace ribll