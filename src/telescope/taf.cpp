#include "include/telescope/taf.h"

#include <TMath.h>

#include "include/event/dssd_event.h"
#include "include/event/csi_event.h"
#include "include/event/telescope_event.h"

namespace ribll {

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


}		// namespace ribll