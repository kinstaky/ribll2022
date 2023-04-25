#include "include/detector/t0d3.h"
#include "include/event/t0_event.h"
#include "include/event/particle_type_event.h"

namespace ribll {

// center of t0d3, in mm
const ROOT::Math::XYZVector t0d3_center{0.0, 0.0, 123.52};
// x range of t0d3, in mm
const std::pair<double, double> t0d3_x_range{-32.0, 32.0};
// y range of t0d3, in mm
const std::pair<double, double> t0d3_y_range{-32.0, 32.0};


T0d3::T0d3(unsigned int run, const std::string &tag)
: Dssd(run, "t0d3", tag) {

	center_ = t0d3_center;
	x_range_ = t0d3_x_range;
	y_range_ = t0d3_y_range;
}


int T0d3::NormalizeSides(TChain *chain, int iteration) {
	if (SideNormalize(chain, 0, 20, iteration)) {
		std::cerr << "Error: Normalize first side failed.\n";
		return -1;
	}
	if (SideNormalize(chain, 1, 13, iteration)) {
		std::cerr << "Error: Normalize second side failed.\n";
		return -1;
	}
	return 0;
}


int T0d3::NormalizeFilter(int iteration) {
	if (iteration == 1) {
		// telescope file name
		TString tele_file_name;
		tele_file_name.Form(
			"%s%st0-telescope-%s%04u.root",
			kGenerateDataPath,
			kTelescopeDir,
			tag_.empty() ? "" : (tag_ + "-").c_str(),
			run_
		);
		// telescope file
		TFile tele_file(tele_file_name, "read");
		// telescope tree
		TTree *ipt = (TTree*)tele_file.Get("tree");
		// telescope event
		T0Event event;
		// setup input branches
		event.SetupInput(ipt);

		// output file name
		TString filter_file_name;
		filter_file_name.Form(
			"%s%st0d3-filter-%s%04u-%d.root",
			kGenerateDataPath,
			kNormalizeDir,
			tag_.empty() ? "" : (tag_ + "-").c_str(),
			run_,
			iteration
		);
		// filter file
		TFile filter_file(filter_file_name, "recreate");
		// filter tree
		TTree opt("tree", "filter tree");
		// filter flag
		unsigned short flag;
		// setup output branch
		opt.Branch("flag", &flag, "f/s");

		// get cut
		std::unique_ptr<TCutG> cut = ReadCut("t0-d2d3-norm-1");
		if (!cut) {
			std::cerr << "Error: Read cut from file failed.\n";
			return -1;
		}

		// filter event, for statistics
		long long filter_num = 0;

		// total number of entries
		long long entries = ipt->GetEntries();
		// 1/100 of entries
		long long entry100 = entries / 100 + 1;
		// show start
		printf("Filtering T0D3 events in iteration %d   0%%", iteration);
		fflush(stdout);
		for (long long entry = 0; entry < entries; ++entry) {
			// show process
			if (entry % entry100 == 0) {
				printf("\b\b\b\b%3lld%%", entry / entry100);
				fflush(stdout);
			}
			// get event
			ipt->GetEntry(entry);
			// initialize flag
			flag = 0;
			if (
				event.num == 1
				&& (event.flag[0] & 0x6) == 0x6
				&& cut->IsInside(event.energy[0][2], event.energy[0][1])
			) {
				flag = 0x1;
				++filter_num;
			}
			opt.Fill();
		}
		// show finish
		printf("\b\b\b\b100%%\n");

		// save tree
		opt.Write();
		// close files
		filter_file.Close();
		tele_file.Close();

		// show statistics
		std::cout << "Filter rate " << filter_num << " / " << entries
			<< "  " << double(filter_num) / double(entries) << "\n";

	} else if (iteration == 2) {
		// telescope file name
		TString tele_file_name;
		tele_file_name.Form(
			"%s%st0-telescope-%s%04u.root",
			kGenerateDataPath,
			kTelescopeDir,
			tag_.empty() ? "" : (tag_ + "-").c_str(),
			run_
		);
		// telescope file
		TFile tele_file(tele_file_name, "read");
		// telescope tree
		TTree *ipt = (TTree*)tele_file.Get("tree");
		// particle type file name
		TString type_file_name;
		type_file_name.Form(
			"%s%st0-particle-type-%s%04u.root",
			kGenerateDataPath,
			kParticleIdentifyDir,
			tag_.empty() ? "" : (tag_ + "-").c_str(),
			run_
		);
		ipt->AddFriend("type=tree", type_file_name);
		// telescope event
		T0Event event;
		// particle type event
		ParticleTypeEvent type;
		// setup input branches
		event.SetupInput(ipt);
		type.SetupInput(ipt, "type.");

		// output file name
		TString filter_file_name;
		filter_file_name.Form(
			"%s%st0d3-filter-%s%04u-%d.root",
			kGenerateDataPath,
			kNormalizeDir,
			tag_.empty() ? "" : (tag_ + "-").c_str(),
			run_,
			iteration
		);
		// filter file
		TFile filter_file(filter_file_name, "recreate");
		// filter tree
		TTree opt("tree", "filter tree");
		// filter flag
		unsigned short flag;
		// setup output branch
		opt.Branch("flag", &flag, "f/s");

		// filter event, for statistics
		long long filter_num = 0;

		// total number of entries
		long long entries = ipt->GetEntries();
		// 1/100 of entries
		long long entry100 = entries / 100 + 1;
		// show start
		printf("Filtering T0D3 events in iteration %d   0%%", iteration);
		fflush(stdout);
		for (long long entry = 0; entry < entries; ++entry) {
			// show process
			if (entry % entry100 == 0) {
				printf("\b\b\b\b%3lld%%", entry / entry100);
				fflush(stdout);
			}
			// get event
			ipt->GetEntry(entry);
			// initialize flag
			flag = 0;
			if (
				event.num == 1
				&& type.charge[0] > 0
				&& type.mass[0] > 0
			) {
				flag = 0x1;
				++filter_num;
			}
			opt.Fill();
		}
		// show finish
		printf("\b\b\b\b100%%\n");

		// save tree
		opt.Write();
		// close files
		filter_file.Close();
		tele_file.Close();

		// show statistics
		std::cout << "Filter rate " << filter_num << " / " << entries
			<< "  " << double(filter_num) / double(entries) << "\n";
	} else {
		std::cerr << "Error: T0D3 normalize filter in iteration "
			<< iteration << " is not implemented yet.\n";
		return -1;
	}
	return 0;
}


}		// namespace ribll
