#include "include/telescope/taf.h"

#include <fstream>

#include <TF1.h>
#include <TFitResult.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>

#include "include/event/dssd_event.h"
#include "include/event/csi_event.h"
#include "include/event/telescope_event.h"

namespace ribll {

// run used for alpha calibration
const unsigned int alpha_calibration_run[6] = {816, 0, 0, 0, 0, 0};
// energy of alpha source, in MeV
const double alpha_energy[3] = {5.157, 5.486, 5.805};
// histogram range of alpha source peaks
const double alpha_hist_range[6][2] = {
	{1700, 2200},
	{0, 60000},
	{0, 60000},
	{0, 60000},
	{0, 60000},
	{0, 60000}
};
// fit range of alpha source peaks
const double alpha_fit_range[6][6] = {
	{1850.0, 1880.0, 1970.0, 2000.0, 2090.0, 2110.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
};


Taf::Taf(unsigned int run, unsigned int index, const std::string &tag)
: Telescope(run, "taf"+std::to_string(index), tag)
, index_(index) {
}


int Taf::Track() {
	std::vector<TString> file_names;
	// add tafd
	file_names.push_back(
		TString::Format(
			"%s%s%s-merge-%s%04u.root",
			kGenerateDataPath,
			kMergeDir,
			(std::string("tafd")+name_[3]).c_str(),
			tag_.empty() ? "" : (tag_+"-").c_str(),
			run_
		)
	);
	// add tafcsi
	file_names.push_back(
		TString::Format(
			"%s%s%s-fundamental-%s%04u.root",
			kGenerateDataPath,
			kFundamentalDir,
			"tafcsi",
			tag_.empty() ? "" : (tag_+"-").c_str(),
			run_
		)
	);

	// first input file
	TFile *ipf = new TFile(file_names[0], "read");
	// first input tree
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< file_names[0] << " failed.\n";
		return -1;
	}
	// add second tree
	if (!ipt->AddFriend("csi=tree", file_names[1])) {
		std::cerr << "Error: Get tree from "
			<< file_names[1] << " falied.\n";
		return -1;
	}

	// input tafd event
	DssdMergeEvent tafd;
	// input tafcsi event
	CircularCsiFundamentalEvent tafcsi;
	// setup input branches
	tafd.SetupInput(ipt);
	tafcsi.SetupInput(ipt, "csi.");


	// output file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-telescope-%s%04u.root",
		kGenerateDataPath,
		kTelescopeDir,
		name_.c_str(),
		tag_.empty() ? "" : (tag_+"-").c_str(),
		run_
	);
	// output file
	TFile *opf = new TFile(output_file_name, "recreate");
	// output tree
	TTree *opt = new TTree("tree", "telescope");
	// output event
	TelescopeEvent tele;
	// setup output branches
	tele.SetupOutput(opt);

	long long total = 0;
	long long match = 0;

	// total number of entries
	long long entries = ipt->GetEntries();
	// 1/100 of entries
	long long entry100 = entries / 100;
	// show start
	printf("Tracking %s telescope   0%%", name_.c_str());
	fflush(stdout);
	for (long long entry = 0; entry < entries; ++entry) {
		// show process
		if (entry % entry100 == 0) {
			printf("\b\b\b\b%3lld%%", entry / entry100);
			fflush(stdout);
		}

		ipt->GetEntry(entry);

		// initialize telescope event
		tele.particle = 0;
		for (size_t i = 0; i < 4; ++i) {
			tele.flag[i] = 0;
			tele.layer[i] = 0;
		}

		if (tafd.hit == 1) {
			++total;

			// record this tafd layer
			++tele.layer[tele.particle];
			tele.flag[tele.particle] |= 1;
			tele.energy[tele.particle][0] = tafd.energy[0];
			tele.radius[tele.particle][0] = tafd.radius[0];
			tele.theta[tele.particle][0] = tafd.theta[0];
			tele.phi[tele.particle][0] = tafd.phi[0];

			// track next layer
			bool conflict = false;
			for (unsigned int i = 0; i < 2; ++i) {
				unsigned int csi_index = index_ * 2 + i;
				if (tafcsi.time[csi_index] < -9e4) continue;
				double max_phi = csi_index < 10
					? 120 - 30 * csi_index
					: 480 - 30 * csi_index;
				max_phi *= TMath::DegToRad();
				double min_phi = csi_index < 10
					? 90 - 30 * csi_index
					: 450 - 30 * csi_index;
				min_phi *= TMath::DegToRad();
				if (tafd.phi[0] > max_phi || tafd.phi[0] < min_phi) continue;

				// check if conflict
				if ((tele.flag[tele.particle] & 2) != 0) {
					conflict = true;
					break;
				}
				// tafcsi goes through the test
				++tele.layer[tele.particle];
				tele.flag[tele.particle] |= 2;
				tele.energy[tele.particle][1] = tafcsi.energy[csi_index];
			}
			if (!conflict) ++tele.particle;
			else {
				tele.flag[tele.particle] = 0;
				tele.layer[tele.particle] = 0;
			}
		}

		if (tele.particle > 0) ++match;

		opt->Fill();
	}
	// show finish
	printf("\b\b\b\b100%%\n");

	// save and close files
	opt->Write();
	opf->Close();
	ipf->Close();

	std::cout << "Match " << name_ << "  "
		<< match << " / " << total
		<< "  " << double(match) / double(total) << "\n";

	return 0;
}


int Taf::Calibrate() {
	if (AlphaCalibrate()) return -1;
	return 0;
}


int Taf::AlphaCalibrate() {
	// calibrate tafd with alpha source
	// input file name
	TString input_file_name;
	input_file_name.Form(
		"%s%s%s-merge-%04u.root",
		kGenerateDataPath,
		kMergeDir,
		(name_.substr(0, 3) + "d" + name_[3]).c_str(),
		alpha_calibration_run[index_]
	);
	// input file
	TFile *ipf = new TFile(input_file_name, "read");
	// input tree
	TTree *ipt = (TTree*)ipf->Get("tree");
	if (!ipt) {
		std::cerr << "Error: Get tree from "
			<< input_file_name << " failed.\n";
		return -1;
	}
	// input data
	DssdMergeEvent event;
	// setup input branches
	event.SetupInput(ipt);

	// output root file name
	TString output_file_name;
	output_file_name.Form(
		"%s%s%s-calibration-%04u.root",
		kGenerateDataPath,
		kCalibrationDir,
		(name_.substr(0, 3) + "d" + name_[3]).c_str(),
		alpha_calibration_run[index_]
	);
	// output file
	TFile *opf = new TFile(output_file_name, "recreate");
	// output energy histogram
	TH1F *he = new TH1F(
		"he", "fit alpha source peaks",
		500, alpha_hist_range[index_][0], alpha_hist_range[index_][1]
	);

	// output parameters txt file
	TString param_file_name;
	param_file_name.Form(
		"%s%s%s-alpha-cali-param.txt",
		kGenerateDataPath,
		kCalibrationDir,
		(name_.substr(0, 3) + "d" + name_[3]).c_str()
	);
	std::ofstream fout(param_file_name.Data());
	if (!fout.good()) {
		std::cerr << "Error: Open file "
			<< param_file_name << " failed.\n";
		return -1;
	}

	// total number of entries
	long long entries = ipt->GetEntries();
	// fill energy to histogram
	for (long long entry = 0; entry < entries; ++entry) {
		ipt->GetEntry(entry);
		if (event.hit != 1) continue;
		he->Fill(event.energy[0]);
	}

	// fit alpha peaks
	TFitResultPtr fit_result[3];
	for (size_t i = 0; i < 3; ++i) {
		TF1 *fpeak = new TF1(
			TString::Format("fpeak%ld", i),
			"gaus",
			alpha_fit_range[index_][2*i],
			alpha_fit_range[index_][2*i+1]
		);
		fpeak->SetParameter(0, 1000);
		fpeak->SetParameter(1, alpha_fit_range[index_][2*i]+10);
		fit_result[i] = he->Fit(fpeak, "RQS+");
	}


	// calibrate
	TGraph *gcali = new TGraph;
	for (size_t i = 0; i < 3; ++i) {
		gcali->AddPoint(fit_result[i]->Parameter(1), alpha_energy[i]);
	}
	TF1 *fcali = new TF1(
		"fcali", "pol1",
		alpha_hist_range[index_][0], alpha_hist_range[index_][1]
	);
	TFitResultPtr cali_result = gcali->Fit(fcali, "QRS+");

	// show calibrate result
	std::cout << cali_result->Parameter(0) << " "
		<< cali_result->Parameter(1) << "\n";

	// output fit results
	fout << cali_result->Parameter(0) << " "
		<< cali_result->Parameter(1) << " "
		<< cali_result->Chi2() << " "
		<< cali_result->Ndf() << " "
		<< cali_result->Chi2() / cali_result->Ndf() << "\n";
	for (size_t i = 0; i < 3; ++i) {
		fout << fit_result[i]->Parameter(0) << " "
			<< fit_result[i]->Parameter(1) << " "
			<< fit_result[i]->Parameter(2) << " "
			<< fit_result[i]->Chi2() << " "
			<< fit_result[i]->Ndf() << " "
			<< fit_result[i]->Chi2() / fit_result[i]->Ndf() << "\n";
	}
	fout.close();

	// save histogram and close files
	he->Write();
	for (size_t i = 0; i < 3; ++i) {
		fit_result[i]->Write(TString::Format("r%ld", i));
	}
	gcali->Write("gcali");
	cali_result->Write("rcali");
	opf->Close();
	// close input file
	ipf->Close();

	return 0;
}



}		// namespace ribll